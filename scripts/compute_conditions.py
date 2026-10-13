#!/usr/bin/env python3
"""
Add satellite and ionosphere conditions to each processed session.

For every tagged session in data/index.json this computes, every STEP_S
seconds through the session:

  Satellites in view   number of healthy satellites above ELEVATION_MASK_DEG
                       at the site, counting only the constellations in that
                       test (GPS only for test G, all four for GREC, ...)
  Ionosphere VTEC      vertical total electron content above the site, in
                       TEC units (1 TECU delays the GPS L1 signal by ~16 cm)

These are the quantities Trimble GNSS Planning charts, computed for the
actual dates of each session from published data, so nothing has to be
copied by hand. The values are stored in the session file next to PDOP, so
the website plots them the same way.

Data sources (public, no account needed):
  satellites  the IGS daily broadcast navigation file (all constellations,
              with health flags), from the BKG data center; posted the
              day after the observations
  ionosphere  global ionosphere maps (IONEX) from UPC Barcelona (15-minute
              rapid maps, about a day after) or ESA (rapid and final)

Until a session's files are published it is marked "pending" and the
website shows the session without these values; the daily run of the
workflow fills them in automatically once the files appear.

Collector times are local clock times (see COLLECTOR_TIME_ZONE); they are
converted to UTC to match the data files.

Needs numpy. Run locally:  python scripts/compute_conditions.py
"""

import argparse
import gzip
import json
import math
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SESSIONS = ROOT / "data" / "sessions"
INDEX = ROOT / "data" / "index.json"
CACHE = ROOT / ".products-cache"          # downloaded files; not committed

# ---------------------------------------------------------------------------
# Settings. Edit here.
# ---------------------------------------------------------------------------

# Approximate antenna location for each site code: latitude and longitude in
# degrees (WGS84) and ellipsoid height in meters. A few kilometers of error
# makes no visible difference to satellite counts or the ionosphere.
SITE_LOCATIONS = {
    "SJER": (37.098, -119.748, 325.0),   # from CA State Plane zone 3 means
    "FRES": (36.813, -119.746, 70.0),    # Fresno State campus (approximate)
    "ELK":  (36.492, -119.788, 42.0),    # from CA State Plane zone 4 means
}

# Time zone of the collector's timestamps. Daylight saving is handled.
COLLECTOR_TIME_ZONE = "America/Los_Angeles"

# Satellites below this elevation angle are not counted (Trimble GNSS
# Planning and most receivers default to 10 degrees).
ELEVATION_MASK_DEG = 10.0

# Spacing of the computed values, in seconds.
STEP_S = 30

# Bump to recompute every session (for example after changing a setting).
CONDITIONS_VERSION = 1

SAT_NAME = "Satellites in view"
ION_NAME = "Ionosphere VTEC"

# Where to download from, best first. {y} year, {d} day of year (3 digits),
# {ymd} yymmdd, {w} GPS week.
NAV_URLS = [
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{d}/BRDC00IGS_R_{y}{d}0000_01D_MN.rnx.gz",
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{d}/BRDC00WRD_R_{y}{d}0000_01D_MN.rnx.gz",
]
ION_URLS = [
    "http://chapman.upc.es/tomion/rapid/{y}/{d}_{ymd}.15min/UPC0OPSRAP_{y}{d}0000_01D_15M_GIM.INX.Z",
    "http://navigation-office.esa.int/products/gnss-products/{w}/ESA0OPSFIN_{y}{d}0000_01D_02H_GIM.INX.gz",
    "http://navigation-office.esa.int/products/gnss-products/{w}/ESA0OPSRAP_{y}{d}0000_01D_01H_GIM.INX.gz",
    "http://chapman.upc.es/tomion/rapid/{y}/{d}_{ymd}.15min/UPC0OPSFIN_{y}{d}0000_01D_02H_GIM.INX.Z",
]
# Precise orbits, used only by --report to check the broadcast positions.
CHECK_SP3_URL = "http://navigation-office.esa.int/products/gnss-products/{w}/ESA0OPSFIN_{y}{d}0000_01D_05M_ORB.SP3.gz"

# ---------------------------------------------------------------------------

GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
GPS_MINUS_UTC = 18                          # seconds, since 2017
BDT_MINUS_GPS = -14                         # BeiDou time = GPS time - 14 s
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563
USER_AGENT = "Caltrans-GNSS-Study site builder (github.com/scpetersonGNSS/Caltrans-GNSS-Study)"


# ---------------------------------------------------------------------------
# Downloading
# ---------------------------------------------------------------------------

_dead_hosts = set()     # hosts that timed out this run; not tried again


def fetch(url, tries=2):
    """Return the bytes at url, or None if it does not exist or can't be reached."""
    host = urllib.parse.urlsplit(url).netloc
    if host in _dead_hosts:
        return None
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except urllib.error.HTTPError as err:
            if err.code in (403, 404, 410):
                return None
            last = err
        except (urllib.error.URLError, TimeoutError, OSError) as err:
            last = err
        time.sleep(2)
    print(f"  could not reach {host} ({last}); skipping it for the rest of this run", flush=True)
    _dead_hosts.add(host)
    return None


def fmt(template, day):
    week = (datetime(day.year, day.month, day.day, tzinfo=timezone.utc) - GPS_EPOCH).days // 7
    return template.format(y=f"{day.year:04d}", d=f"{day.timetuple().tm_yday:03d}",
                           ymd=day.strftime("%y%m%d"), w=f"{week:04d}")


def load_text(url):
    """Download (or reuse) a compressed product; return its text, or None."""
    CACHE.mkdir(exist_ok=True)
    local = CACHE / Path(urllib.parse.urlsplit(url).path).name
    raw = local.read_bytes() if local.exists() else fetch(url)
    if raw is None:
        return None
    try:
        if url.endswith(".gz"):
            text = gzip.decompress(raw)
        elif url.endswith(".Z"):                      # unix compress; gzip can read it
            text = subprocess.run(["gzip", "-dc"], input=raw, capture_output=True, check=True).stdout
        else:
            text = raw
    except (OSError, EOFError, subprocess.CalledProcessError):
        local.unlink(missing_ok=True)
        return None
    if not local.exists():
        local.write_bytes(raw)
    return text.decode("ascii", "replace")


# ---------------------------------------------------------------------------
# Broadcast navigation (RINEX 3 mixed) -> satellite positions
# ---------------------------------------------------------------------------

def _num(s):
    s = s.strip().replace("D", "E").replace("d", "e")
    return float(s) if s else math.nan


# Lines per record after the first, by system.
NAV_LINES = {"G": 7, "E": 7, "C": 7, "J": 7, "I": 7, "R": 3, "S": 3}


def parse_nav(text):
    """{'G': [record, ...], 'R': [...], ...}. Each record is a dict."""
    out = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines) and "END OF HEADER" not in lines[i]:
        i += 1
    i += 1
    while i < len(lines):
        line = lines[i]
        sysc = line[:1]
        if sysc not in NAV_LINES or len(line) < 23:
            i += 1
            continue
        n = NAV_LINES[sysc]
        block = lines[i:i + n + 1]
        i += n + 1
        try:
            prn = f"{sysc}{int(line[1:3]):02d}"
            y, mo, d, h, mi, s = (int(v) for v in line[3:23].split())
            vals = [_num(line[23 + 19 * k:42 + 19 * k]) for k in range(3)]
            for b in block[1:]:
                vals += [_num(b[4 + 19 * k:23 + 19 * k]) for k in range(4)]
        except (ValueError, IndexError):
            continue
        epoch = datetime(y, mo, d, h, mi, s, tzinfo=timezone.utc)
        out.setdefault(sysc, []).append({"prn": prn, "epoch": epoch, "v": vals})
    return out


# Kepler constants per system: (GM, Earth rotation rate)
KEPLER = {"G": (3.986005e14, 7.2921151467e-5), "J": (3.986005e14, 7.2921151467e-5),
          "E": (3.986004418e14, 7.2921151467e-5), "C": (3.986004418e14, 7.292115e-5)}
BDS_GEO = {1, 2, 3, 4, 5, 59, 60, 61, 62, 63}


def kepler_positions(rec, t_sys):
    """ECEF positions (m) of a GPS / Galileo / BeiDou / QZSS satellite at
    system times t_sys (seconds of the same week numbering as Toe)."""
    v = rec["v"]
    sysc = rec["prn"][0]
    gm, we = KEPLER[sysc]
    crs, dn, m0 = v[4], v[5], v[6]
    cuc, e, cus, sqa = v[7], v[8], v[9], v[10]
    toe, cic, om0, cis = v[11], v[12], v[13], v[14]
    i0, crc, w, omd = v[15], v[16], v[17], v[18]
    idot = v[19]
    a = sqa * sqa
    tk = t_sys - toe
    tk = (tk + 302400) % 604800 - 302400       # week crossover
    n = math.sqrt(gm / a ** 3) + dn
    M = m0 + n * tk
    E = M.copy()
    for _ in range(10):
        E = M + e * np.sin(E)
    nu = np.arctan2(math.sqrt(1 - e * e) * np.sin(E), np.cos(E) - e)
    phi = nu + w
    du = cus * np.sin(2 * phi) + cuc * np.cos(2 * phi)
    dr = crs * np.sin(2 * phi) + crc * np.cos(2 * phi)
    di = cis * np.sin(2 * phi) + cic * np.cos(2 * phi)
    u = phi + du
    r = a * (1 - e * np.cos(E)) + dr
    inc = i0 + di + idot * tk
    xo, yo = r * np.cos(u), r * np.sin(u)
    geo = sysc == "C" and int(rec["prn"][1:]) in BDS_GEO
    if geo:
        om = om0 + omd * tk - we * toe
    else:
        om = om0 + (omd - we) * tk - we * toe
    x = xo * np.cos(om) - yo * np.cos(inc) * np.sin(om)
    y = xo * np.sin(om) + yo * np.cos(inc) * np.cos(om)
    z = yo * np.sin(inc)
    if geo:                                    # BeiDou GEO: rotate from the -5° frame
        f, g = we * tk, math.radians(-5)
        xg, yg, zg = x, y * math.cos(g) + z * math.sin(g), -y * math.sin(g) + z * math.cos(g)
        x = xg * np.cos(f) + yg * np.sin(f)
        y = -xg * np.sin(f) + yg * np.cos(f)
        z = zg
    return np.column_stack([x, y, z])


GLO_GM, GLO_AE, GLO_J2, GLO_WE = 398600.4418e9, 6378136.0, 1.08262575e-3, 7.292115e-5


def _glo_deriv(s, acc):
    x, y, z, vx, vy, vz = s
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    a = GLO_GM / r ** 3
    b = 1.5 * GLO_J2 * GLO_GM * GLO_AE ** 2 / r ** 5
    zz = 5 * z * z / r2
    return np.array([vx, vy, vz,
                     -a * x - b * x * (1 - zz) + GLO_WE ** 2 * x + 2 * GLO_WE * vy + acc[0],
                     -a * y - b * y * (1 - zz) + GLO_WE ** 2 * y - 2 * GLO_WE * vx + acc[1],
                     -a * z - b * z * (3 - zz) + acc[2]])


def glonass_positions(rec, dt):
    """ECEF positions (m) at dt seconds from the record epoch (|dt| <= ~20 min),
    by RK4 integration of the broadcast state (GLONASS ICD)."""
    v = rec["v"]
    s0 = np.array([v[3], v[7], v[11], v[4], v[8], v[12]]) * 1000.0
    acc = np.array([v[5], v[9], v[13]]) * 1000.0
    step = 30.0
    out = np.empty((len(dt), 3))
    for sign in (1, -1):
        sel = np.where((dt >= 0) if sign > 0 else (dt < 0))[0]
        if not len(sel):
            continue
        order = sel[np.argsort(np.abs(dt[sel]))]
        s, t = s0.copy(), 0.0
        for k in order:
            target = abs(dt[k])
            while t < target - 1e-9:
                hh = min(step, target - t) * sign
                k1 = _glo_deriv(s, acc); k2 = _glo_deriv(s + hh / 2 * k1, acc)
                k3 = _glo_deriv(s + hh / 2 * k2, acc); k4 = _glo_deriv(s + hh * k3, acc)
                s = s + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                t += abs(hh)
            out[k] = s[:3]
    return out


def healthy(rec):
    v, sysc = rec["v"], rec["prn"][0]
    if sysc == "R":
        return v[6] == 0                          # Bn health flag
    if sysc in ("G", "J"):
        return v[24] == 0                         # SV health
    if sysc == "E":
        return v[24] == 0                         # SV health bits (E1B/E5a/E5b)
    if sysc == "C":
        return v[24] == 0                         # SatH1
    return True


def broadcast_positions(nav, letters, tq):
    """{prn: positions Nx3 (NaN where no valid ephemeris)} at UTC posix times tq."""
    gps_sow = (tq + GPS_MINUS_UTC - GPS_EPOCH.timestamp()) % 604800
    out = {}
    for sysc in letters:
        recs = [r for r in nav.get(sysc, []) if healthy(r) and np.all(np.isfinite(r["v"][:19 if sysc != "R" else 15]))]
        by_prn = {}
        for r in recs:
            by_prn.setdefault(r["prn"], []).append(r)
        for prn, rs in by_prn.items():
            if sysc == "R":
                ref = np.array([r["epoch"].timestamp() for r in rs])                       # UTC
                limit = 20 * 60
            else:
                # Toe as UTC posix (approx.; only used to pick the nearest record).
                ref = np.array([r["epoch"].timestamp() - (GPS_MINUS_UTC if sysc != "C" else GPS_MINUS_UTC + BDT_MINUS_GPS) for r in rs])
                limit = 4 * 3600 if sysc in ("G", "J") else 3 * 3600
            pick = np.abs(tq[:, None] - ref[None, :]).argmin(axis=1)
            near = np.abs(tq - ref[pick]) <= limit
            pos = np.full((len(tq), 3), np.nan)
            for j in np.unique(pick[near]):
                idx = np.where(near & (pick == j))[0]
                if sysc == "R":
                    pos[idx] = glonass_positions(rs[j], tq[idx] - ref[j])
                else:
                    t_sys = gps_sow[idx] + (BDT_MINUS_GPS if sysc == "C" else 0)
                    pos[idx] = kepler_positions(rs[j], t_sys)
            out[prn] = pos
    return out


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def station_ecef(lat, lon, h):
    la, lo = math.radians(lat), math.radians(lon)
    e2 = WGS84_F * (2 - WGS84_F)
    N = WGS84_A / math.sqrt(1 - e2 * math.sin(la) ** 2)
    r = np.array([(N + h) * math.cos(la) * math.cos(lo), (N + h) * math.cos(la) * math.sin(lo), (N * (1 - e2) + h) * math.sin(la)])
    east = np.array([-math.sin(lo), math.cos(lo), 0.0])
    north = np.array([-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la)])
    up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
    return r, np.vstack([east, north, up])


def geometry(pos, site, want_pdop=False):
    """Satellites above the mask at each time, and optionally PDOP (one
    receiver clock per constellation, as multi-GNSS receivers solve)."""
    r0, enu = station_ecef(*site)
    n = len(next(iter(pos.values()))) if pos else 0
    count = np.zeros(n, dtype=int)
    rows = []
    systems = sorted({p[0] for p in pos})
    for prn, p in pos.items():
        los = p - r0
        rng = np.linalg.norm(los, axis=1)
        u = (los @ enu.T) / rng[:, None]
        el = np.degrees(np.arcsin(np.clip(np.nan_to_num(u[:, 2], nan=-1), -1, 1)))
        vis = np.isfinite(rng) & (el >= ELEVATION_MASK_DEG)
        count += vis
        rows.append((vis, u, systems.index(prn[0])))
    if not want_pdop:
        return count, None
    pdop = np.full(n, np.nan)
    for k in range(n):
        G = [list(-u[k]) + [1.0 if s == j else 0.0 for j in range(len(systems))] for vis, u, s in rows if vis[k]]
        if not G:
            continue
        G = np.array(G)
        G = np.hstack([G[:, :3], G[:, 3:][:, G[:, 3:].any(axis=0)]])
        if len(G) < G.shape[1]:
            continue
        try:
            Q = np.linalg.inv(G.T @ G)
        except np.linalg.LinAlgError:
            continue
        pdop[k] = math.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])
    return count, pdop


# ---------------------------------------------------------------------------
# Ionosphere (IONEX)
# ---------------------------------------------------------------------------

def parse_ionex(text):
    """List of (time UTC posix, lats, lons, grid[lat, lon] in TECU)."""
    head, _, body = text.partition("END OF HEADER")
    exponent = -1
    lat1 = dlat = lon1 = lon2 = dlon = None
    for line in head.splitlines():
        label = line[60:].strip()
        if label == "EXPONENT":
            exponent = int(line[:6])
        elif label == "LAT1 / LAT2 / DLAT":
            lat1, lat2, dlat = map(float, line[:60].split()[:3])
        elif label == "LON1 / LON2 / DLON":
            lon1, lon2, dlon = map(float, line[:60].split()[:3])
    if dlat is None or dlon is None:
        return []
    lats = np.arange(lat1, lat2 + dlat / 2, dlat)
    lons = np.arange(lon1, lon2 + dlon / 2, dlon)
    nlon = len(lons)
    maps, grid, epoch, scale = [], None, None, 10.0 ** exponent
    lines = body.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        label = line[60:].strip()
        if label in ("START OF RMS MAP", "START OF HEIGHT MAP"):
            break
        if label == "START OF TEC MAP":
            grid = np.full((len(lats), nlon), np.nan); epoch = None
        elif label == "EPOCH OF CURRENT MAP":
            y, mo, d, h, mi, s = map(int, line[:36].split())
            epoch = (datetime(y, mo, d, tzinfo=timezone.utc) + timedelta(hours=h, minutes=mi, seconds=s)).timestamp()
        elif label == "EXPONENT":
            scale = 10.0 ** int(line[:6])
        elif label == "LAT/LON1/LON2/DLON/H" and grid is not None:
            lat = float(line[2:8])
            vals = []
            while len(vals) < nlon and i + 1 < len(lines):
                i += 1
                row = lines[i].rstrip()
                vals += [int(row[k:k + 5]) for k in range(0, len(row), 5) if row[k:k + 5].strip()]
            j = int(round((lat - lat1) / dlat))
            v = np.array(vals[:nlon], dtype=float)
            v[v == 9999] = np.nan
            if 0 <= j < len(lats) and len(v) == nlon:
                grid[j] = v * scale
        elif label == "END OF TEC MAP" and grid is not None and epoch is not None:
            maps.append((epoch, lats, lons, grid))
            grid = None
        i += 1
    return maps


def vtec_at(maps, tq, lat, lon):
    maps = sorted({m[0]: m for m in maps}.values(), key=lambda m: m[0])
    times = np.array([m[0] for m in maps])

    def at_map(m):
        _, lats, lons, g = m
        fy = (lat - lats[0]) / (lats[1] - lats[0]); fx = (lon - lons[0]) / (lons[1] - lons[0])
        j = min(max(int(math.floor(fy)), 0), len(lats) - 2); i = min(max(int(math.floor(fx)), 0), len(lons) - 2)
        wy, wx = fy - j, fx - i
        return ((1 - wy) * ((1 - wx) * g[j, i] + wx * g[j, i + 1]) + wy * ((1 - wx) * g[j + 1, i] + wx * g[j + 1, i + 1]))

    vals = np.array([at_map(m) for m in maps])
    out = np.interp(tq, times, vals, left=np.nan, right=np.nan)
    gi = np.clip(np.searchsorted(times, tq), 1, len(times) - 1)
    out[(times[gi] - times[gi - 1]) > 3 * 3600 + 1] = np.nan
    return out


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def utc_days(t0, t1, margin_s=0):
    d = datetime.fromtimestamp(t0 - margin_s, timezone.utc).date()
    last = datetime.fromtimestamp(t1 + margin_s, timezone.utc).date()
    days = []
    while d <= last:
        days.append(d); d += timedelta(days=1)
    return days


def first_available(urls, day, parse):
    for template in urls:
        url = fmt(template, day)
        text = load_text(url)
        if text:
            parsed = parse(text)
            if parsed:
                return parsed, Path(urllib.parse.urlsplit(url).path).name
    return None, None


def get_nav(days):
    nav, used = {}, []
    for day in days:
        got, name = first_available(NAV_URLS, day, parse_nav)
        if not got:
            return None, f"no broadcast navigation file yet for {day}"
        for k, v in got.items():
            nav.setdefault(k, []).extend(v)
        used.append(name)
    return (nav, used), None


def get_ionosphere(days):
    maps, used = [], []
    for day in days:
        got, name = first_available(ION_URLS, day, parse_ionex)
        if not got:
            return None, f"no ionosphere map yet for {day}"
        maps += got; used.append(name)
    return (maps, used), None


def session_grid(data, utc_offset_h=None):
    """UTC posix times of the computation grid, and the grid in s from start."""
    start = datetime.fromisoformat(data["start"])
    tz = ZoneInfo(COLLECTOR_TIME_ZONE) if utc_offset_h is None else timezone(timedelta(hours=utc_offset_h))
    start_utc = start.replace(tzinfo=tz).timestamp()
    end = data["t"][-1]
    grid = np.arange(0, end + STEP_S / 2, STEP_S, dtype=float)
    return start_utc + grid, grid


def as_extra(values, dec, unit, note):
    scale = 10 ** dec
    return {"unit": unit, "dec": dec, "scale": scale, "step_s": STEP_S, "computed": True, "note": note,
            "v": [None if not np.isfinite(v) else int(round(v * scale)) for v in values]}


def compute(summary, data):
    site = SITE_LOCATIONS[summary["site"]]
    letters = summary["test"]
    tq, _ = session_grid(data)
    nav, why = get_nav(utc_days(tq[0], tq[-1], 4 * 3600))
    if not nav:
        return None, why
    ion, why = get_ionosphere(utc_days(tq[0], tq[-1]))
    if not ion:
        return None, why
    (nav, nav_files), (maps, ion_files) = nav, ion
    pos = broadcast_positions(nav, letters, tq)
    count, _ = geometry(pos, site)
    vtec = vtec_at(maps, tq, site[0], site[1])
    cons = " + ".join(summary.get("constellations", [])) or letters
    extras = {
        SAT_NAME: as_extra(count.astype(float), 0, "",
                           f"Satellites in view are computed, not logged: healthy {cons} satellites above "
                           f"{ELEVATION_MASK_DEG:g}° at the site, from the broadcast orbits ({', '.join(nav_files)})."),
        ION_NAME: as_extra(vtec, 1, "TECU",
                           f"Ionosphere VTEC is computed, not logged: vertical total electron content above the site "
                           f"from global ionosphere maps ({', '.join(ion_files)}). 1 TECU ≈ 16 cm of GPS L1 delay."),
    }
    info = {"status": "done", "version": CONDITIONS_VERSION, "mask_deg": ELEVATION_MASK_DEG, "step_s": STEP_S,
            "time_zone": COLLECTOR_TIME_ZONE, "navigation": nav_files, "ionosphere": ion_files,
            "computed": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    return (extras, info, nav), None


# ---------------------------------------------------------------------------
# Checks (--report)
# ---------------------------------------------------------------------------

def parse_sp3(text):
    sats, t, gps = {}, None, True
    for line in text.splitlines():
        if line.startswith("%c") and len(line) > 13:
            gps = line[9:12].strip() in ("GPS", "ccc", "")
        elif line.startswith("*"):
            p = line[1:].split()
            t = datetime(*map(int, p[:5]), tzinfo=timezone.utc).timestamp() + float(p[5]) - (GPS_MINUS_UTC if gps else 0)
        elif line.startswith("P") and t is not None and len(line) >= 46:
            try:
                xyz = tuple(float(line[4 + 14 * k:18 + 14 * k]) * 1000 for k in range(3))
            except ValueError:
                continue
            if xyz != (0.0, 0.0, 0.0):
                sats.setdefault(line[1:4].replace(" ", "0"), []).append((t,) + xyz)
    return {p: np.array(sorted(v)) for p, v in sats.items()}


def recorded_on_grid(data, grid):
    x = (data.get("extra") or {}).get("PDOP")
    if not x:
        return None
    t = np.array(data["t"], dtype=float)
    v = np.array([np.nan if a is None else a / x["scale"] for a in x["v"]])
    i = np.clip(np.searchsorted(t, grid), 1, len(t) - 1)
    i = np.where(np.abs(t[i - 1] - grid) < np.abs(t[i] - grid), i - 1, i)
    out = v[i]
    out[np.abs(t[i] - grid) > 5] = np.nan
    return out


def report_checks(summary, data, nav, lines):
    site, letters = SITE_LOCATIONS[summary["site"]], summary["test"]
    # 1. Broadcast vs precise orbit positions, at the precise orbit's epochs.
    day = datetime.fromisoformat(data["start"]).date() + timedelta(days=1)
    sp3 = load_text(fmt(CHECK_SP3_URL, day))
    if sp3:
        prec = parse_sp3(sp3)
        errs = {}
        for prn, a in prec.items():
            if prn[0] not in letters:
                continue
            tq = a[::12, 0]
            p = broadcast_positions(nav, prn[0], tq).get(prn)
            if p is None:
                continue
            d = np.linalg.norm(p - a[::12, 1:4], axis=1)
            d = d[np.isfinite(d)]
            if len(d):
                errs.setdefault(prn[0], []).append(float(np.median(d)))
        lines.append("- broadcast vs ESA precise orbit, median error per satellite (m): " +
                     "; ".join(f"{k}: median {np.median(v):.1f}, worst {max(v):.1f} ({len(v)} sats)" for k, v in sorted(errs.items())))
    else:
        lines.append("- precise orbit check skipped (ESA file not found)")
    # 2. Time zone: correlation of computed vs recorded PDOP for several offsets.
    if "PDOP" in (data.get("extra") or {}):
        res = []
        for off in (-9, -8, -7, -6, 0):
            tq, grid = session_grid(data, off)
            sub = slice(None, None, 4)                    # every 2 minutes is plenty
            _, pdop = geometry(broadcast_positions(nav, letters, tq[sub]), site, want_pdop=True)
            rec = recorded_on_grid(data, grid[sub])
            ok = np.isfinite(pdop) & np.isfinite(rec)
            r = float(np.corrcoef(pdop[ok], rec[ok])[0, 1]) if ok.sum() > 20 else float("nan")
            res.append(f"UTC{off:+d}h r={r:.3f}")
            if off == -7:
                lines.append(f"- PDOP at UTC-7: recorded mean {np.nanmean(rec):.2f}, computed mean {np.nanmean(pdop):.2f}")
        lines.append("- recorded vs computed PDOP correlation by collector time offset: " + ", ".join(res))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", help="comma-separated session ids to compute")
    ap.add_argument("--force", action="store_true", help="recompute even if already done")
    ap.add_argument("--report", help="also run checks and write them (Markdown) to this file")
    args = ap.parse_args()

    index = json.loads(INDEX.read_text())
    only = set(args.only.split(",")) if args.only else None
    report = ["# Conditions check", ""]
    changed = False
    for s in index.get("sessions", []):
        if only and s["id"] not in only:
            continue
        if s.get("site") not in SITE_LOCATIONS or not s.get("test"):
            continue
        path = SESSIONS / f"{s['id']}.json"
        data = json.loads(path.read_text())
        extra = data.setdefault("extra", {})
        cond = s.get("conditions") or {}
        current = (SAT_NAME in extra and ION_NAME in extra and cond.get("status") == "done"
                   and cond.get("version") == CONDITIONS_VERSION)
        if current and not args.force and not args.report:
            continue
        print(f"{s['id']} ({s['site']} / {s['test']}):", flush=True)
        result, why = compute(s, data)
        if not result:
            print(f"  pending: {why}", flush=True)
            if cond.get("status") != "done":
                s["conditions"] = {"status": "pending", "reason": why, "version": CONDITIONS_VERSION}
                changed |= cond != s["conditions"]
            report += [f"## {s['id']} {s['site']} {s['test']}", f"- pending: {why}", ""]
            continue
        extras, info, nav = result
        sats = [v for v in extras[SAT_NAME]["v"] if v is not None]
        ions = [v / 10 for v in extras[ION_NAME]["v"] if v is not None]
        print(f"  satellites {min(sats)}–{max(sats)}, VTEC {min(ions):.1f}–{max(ions):.1f} TECU "
              f"from {', '.join(info['navigation'])} and {', '.join(info['ionosphere'])}", flush=True)
        if args.report:
            sv, iv = extras[SAT_NAME]["v"], extras[ION_NAME]["v"]
            lines = [f"## {s['id']} {s['site']} {s['test']} (start {data['start']} local)",
                     f"- files: {', '.join(info['navigation'])}; {', '.join(info['ionosphere'])}",
                     f"- satellites: min {min(sats)}, max {max(sats)}, mean {sum(sats) / len(sats):.1f}",
                     f"- VTEC: min {min(ions):.1f}, max {max(ions):.1f} TECU",
                     "- every 2 h from start (sats / TECU): " + ", ".join(
                         f"{sv[k]}/{'-' if iv[k] is None else iv[k] / 10}" for k in range(0, len(sv), 240))]
            report_checks(s, data, nav, lines)
            report += lines + [""]
        extra.update(extras)
        path.write_text(json.dumps(data, separators=(",", ":")))
        s["conditions"] = info
        changed = True

    if args.report:
        Path(args.report).write_text("\n".join(report) + "\n")
    if changed:
        INDEX.write_text(json.dumps(index, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
