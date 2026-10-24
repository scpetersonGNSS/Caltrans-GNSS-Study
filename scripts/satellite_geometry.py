#!/usr/bin/env python3
"""
Predict satellite visibility and DOP for each processed session, in the style
of Trimble's GNSS Planning charts (number of satellites, DOP over time).

For every session in data/index.json this script:
  1. finds where the receiver was (session mean Northing / Easting, converted
     from the site's state plane zone to latitude / longitude),
  2. downloads the IGS broadcast orbits (RINEX 3 mixed navigation file) for
     each UTC day the session covers, caching them in data/nav/,
  3. computes, every STEP_S seconds, every satellite's elevation and azimuth,
     counts the healthy satellites above ELEV_MASK_DEG for each constellation,
     and computes PDOP / HDOP / VDOP for the constellations used in that test,
  4. writes data/geometry/<session-id>.json for the website.

Sessions that already have an up-to-date geometry file are skipped. Sessions
without a site tag are skipped (their location is not known).

Uses only the Python standard library.

Run locally:   python scripts/satellite_geometry.py
Rebuild all:   python scripts/satellite_geometry.py --all
"""

import gzip
import json
import math
import sys
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parent.parent
INDEX = ROOT / "data" / "index.json"
GEOMETRY = ROOT / "data" / "geometry"
NAV = ROOT / "data" / "nav"

# ---------------------------------------------------------------------------
# Settings. Edit here.
# ---------------------------------------------------------------------------

# Satellites below this elevation are not counted (Trimble GNSS Planning's
# default cutoff is 10 degrees).
ELEV_MASK_DEG = 10.0

# Time between predicted epochs, in seconds.
STEP_S = 60

# Time zone of the timestamps in the uploaded files.
DATA_TIME_ZONE = "America/Los_Angeles"

# Coordinate system of each site's Northing / Easting (Lambert conformal conic,
# meters). California zone 3 = NAD83 / California zone 3 (EPSG:26943).
CA_ZONE_3 = {"lat1": 37 + 4 / 60, "lat2": 38 + 26 / 60, "lat0": 36.5,
             "lon0": -120.5, "fe": 2000000.0, "fn": 500000.0}
CA_ZONE_4 = {"lat1": 36.0, "lat2": 37.25, "lat0": 35 + 20 / 60,
             "lon0": -119.0, "fe": 2000000.0, "fn": 500000.0}
SITE_PROJECTION = {"SJER": CA_ZONE_3, "ELK": CA_ZONE_4, "FRES": CA_ZONE_4}

# Where to download broadcast orbits from, tried in order. {y} year, {j} day
# of year (3 digits). All are public and need no login.
NAV_SOURCES = [
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{j}/BRDC00WRD_S_{y}{j}0000_01D_MN.rnx.gz",
    "https://igs.bkg.bund.de/root_ftp/MGEX/BRDC/{y}/{j}/BRDM00DLR_S_{y}{j}0000_01D_MN.rnx.gz",
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{j}/BRDC00IGS_R_{y}{j}0000_01D_MN.rnx.gz",
    "https://gssc.esa.int/igs/data/daily/{y}/{j}/BRDC00WRD_S_{y}{j}0000_01D_MN.rnx.gz",
]

# Constellation letters, in chart order.
SYSTEMS = "GREC"

# Bump when the geometry file layout or the computation changes; older files
# are rebuilt automatically.
FORMAT_VERSION = 1

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
LEAP_SECONDS = 18                 # GPS - UTC since 2017
BDT_WEEK_OFFSET = 1356            # BeiDou week 0 = GPS week 1356
BDT_GPST_OFFSET = 14              # BDT = GPST - 14 s

WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563
GRS80_F = 1 / 298.257222101

MU = {"G": 3.9860050e14, "E": 3.986004418e14, "C": 3.986004418e14}
OMEGA_E = {"G": 7.2921151467e-5, "E": 7.2921151467e-5, "C": 7.292115e-5}
MAX_DTOE = {"G": 7200, "E": 14400, "C": 21600, "R": 1800}

GLO_MU = 3.9860044e14
GLO_J2 = 1.0826257e-3
GLO_RE = 6378136.0
GLO_OMEGA = 7.292115e-5

SIN_5 = -0.0871557427476582       # sin(-5 deg), BeiDou GEO rotation
COS_5 = 0.9961946980917456


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------

def lcc_inverse(northing, easting, p, f=GRS80_F, a=WGS84_A):
    """Lambert conformal conic (two standard parallels) -> lat, lon degrees."""
    e2 = f * (2 - f)
    e = math.sqrt(e2)

    def m(phi):
        return math.cos(phi) / math.sqrt(1 - e2 * math.sin(phi) ** 2)

    def t(phi):
        s = e * math.sin(phi)
        return math.tan(math.pi / 4 - phi / 2) / ((1 - s) / (1 + s)) ** (e / 2)

    p1, p2, p0 = (math.radians(p[k]) for k in ("lat1", "lat2", "lat0"))
    n = (math.log(m(p1)) - math.log(m(p2))) / (math.log(t(p1)) - math.log(t(p2)))
    F = m(p1) / (n * t(p1) ** n)
    rho0 = a * F * t(p0) ** n

    x = easting - p["fe"]
    y = rho0 - (northing - p["fn"])
    rho = math.copysign(math.hypot(x, y), n)
    tt = (rho / (a * F)) ** (1 / n)
    theta = math.atan2(x, y)
    lon = theta / n + math.radians(p["lon0"])
    phi = math.pi / 2 - 2 * math.atan(tt)
    for _ in range(10):
        s = e * math.sin(phi)
        phi = math.pi / 2 - 2 * math.atan(tt * ((1 - s) / (1 + s)) ** (e / 2))
    return math.degrees(phi), math.degrees(lon)


def geodetic_to_ecef(lat, lon, h):
    e2 = WGS84_F * (2 - WGS84_F)
    la, lo = math.radians(lat), math.radians(lon)
    N = WGS84_A / math.sqrt(1 - e2 * math.sin(la) ** 2)
    return ((N + h) * math.cos(la) * math.cos(lo),
            (N + h) * math.cos(la) * math.sin(lo),
            (N * (1 - e2) + h) * math.sin(la))


def enu_rotation(lat, lon):
    la, lo = math.radians(lat), math.radians(lon)
    sl, cl, so, co = math.sin(la), math.cos(la), math.sin(lo), math.cos(lo)
    return ((-so, co, 0.0),
            (-sl * co, -sl * so, cl),
            (cl * co, cl * so, sl))


# ---------------------------------------------------------------------------
# Time
# ---------------------------------------------------------------------------

def utc_to_gps_seconds(dt):
    return (dt - GPS_EPOCH).total_seconds() + LEAP_SECONDS


def local_to_utc(text, tz):
    return datetime.fromisoformat(text).replace(tzinfo=ZoneInfo(tz)).astimezone(timezone.utc)


# ---------------------------------------------------------------------------
# Broadcast orbits
# ---------------------------------------------------------------------------

def nav_file_for(day):
    """Return the cached navigation file for a UTC date, downloading it if needed."""
    y, j = f"{day.year:04d}", f"{day.timetuple().tm_yday:03d}"
    NAV.mkdir(parents=True, exist_ok=True)
    cached = sorted(NAV.glob(f"*_{y}{j}0000_01D_MN.rnx.gz"))
    if cached:
        return cached[0]
    errors = []
    for template in NAV_SOURCES:
        url = template.format(y=y, j=j)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "caltrans-gnss-study"})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            gzip.decompress(body)[:100]            # make sure it is a real gzip file
        except Exception as err:                   # noqa: BLE001 - try the next source
            errors.append(f"{url}: {err}")
            continue
        path = NAV / url.rsplit("/", 1)[1]
        path.write_bytes(body)
        print(f"Downloaded {path.name} ({len(body) // 1024} kB)")
        return path
    raise RuntimeError(f"No broadcast orbits found for {day}:\n  " + "\n  ".join(errors))


def _num(text):
    text = text.strip().replace("D", "E").replace("d", "e")
    return float(text) if text else 0.0


def parse_nav(path):
    """Read a RINEX 3 navigation file into {sat: [record, ...]}."""
    lines = gzip.decompress(path.read_bytes()).decode("ascii", errors="replace").splitlines()
    i = 0
    while i < len(lines) and "END OF HEADER" not in lines[i]:
        i += 1
    i += 1
    extra_lines = {"G": 7, "E": 7, "C": 7, "J": 7, "I": 7, "R": 3, "S": 3}
    version = 3.0
    for line in lines[:i]:
        if "RINEX VERSION" in line:
            version = float(line[:9])
    if version >= 3.05:
        extra_lines["R"] = 4
    records = {}
    while i < len(lines):
        line = lines[i]
        if not line.strip() or line[0] not in extra_lines or not line[1:3].strip().isdigit():
            i += 1
            continue
        sys_id = line[0]
        count = extra_lines[sys_id]
        block = lines[i:i + 1 + count]
        i += 1 + count
        if sys_id not in SYSTEMS:
            continue
        try:
            sat = f"{sys_id}{int(line[1:3]):02d}"
            ep = [int(v) for v in line[4:23].split()]
            toc = datetime(ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], tzinfo=timezone.utc)
            vals = [_num(line[23 + 19 * k:42 + 19 * k]) for k in range(3)]
            for cont in block[1:]:
                vals += [_num(cont[4 + 19 * k:23 + 19 * k]) for k in range(4)]
        except (ValueError, IndexError):
            continue
        rec = make_record(sys_id, toc, vals)
        if rec:
            records.setdefault(sat, []).append(rec)
    return records


# The RINEX 3 Keplerian layout (values after the epoch line's three clock terms):
#  0 IODE  1 Crs  2 dn   3 M0
#  4 Cuc   5 e    6 Cus  7 sqrtA
#  8 toe   9 Cic 10 OMEGA0 11 Cis
# 12 i0   13 Crc 14 omega 15 OMEGADOT
# 16 IDOT 17 codes/spare 18 week 19 flag/spare
# 20 accuracy 21 health 22 TGD 23 IODC
# vals[] passed to make_record starts with the three clock terms, so add 3.
def make_record(sys_id, toc, v):
    """One broadcast ephemeris record; None if it cannot be used."""
    if sys_id == "R":
        return {"t": utc_to_gps_seconds(toc), "healthy": int(v[3 + 3]) == 0,
                "pos": [v[3] * 1e3, v[7] * 1e3, v[11] * 1e3],
                "vel": [v[4] * 1e3, v[8] * 1e3, v[12] * 1e3],
                "acc": [v[5] * 1e3, v[9] * 1e3, v[13] * 1e3]}
    if len(v) < 3 + 24:
        return None
    k = v[3:]
    sqa, toe, week = k[7], k[8], k[18]
    if sqa <= 0:
        return None
    if sys_id == "C":
        t_toe = (week + BDT_WEEK_OFFSET) * 604800 + toe + BDT_GPST_OFFSET
    else:
        t_toe = week * 604800 + toe
    return {"t": t_toe, "healthy": int(k[21]) == 0,
            "crs": k[1], "dn": k[2], "m0": k[3], "cuc": k[4], "e": k[5], "cus": k[6],
            "sqa": sqa, "toe": toe, "cic": k[9], "om0": k[10], "cis": k[11], "i0": k[12],
            "crc": k[13], "w": k[14], "omd": k[15], "idot": k[16]}


def kepler_position(sys_id, prn, r, t):
    """ECEF position (m) of a GPS / Galileo / BeiDou satellite at GPS time t."""
    mu, we = MU[sys_id], OMEGA_E[sys_id]
    tk = t - r["t"]
    a = r["sqa"] ** 2
    n = math.sqrt(mu / a ** 3) + r["dn"]
    M = r["m0"] + n * tk
    e = r["e"]
    E = M
    for _ in range(30):
        dE = (M - E + e * math.sin(E)) / (1 - e * math.cos(E))
        E += dE
        if abs(dE) < 1e-13:
            break
    nu = math.atan2(math.sqrt(1 - e * e) * math.sin(E), math.cos(E) - e)
    phi = nu + r["w"]
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)
    u = phi + r["cus"] * s2 + r["cuc"] * c2
    rad = a * (1 - e * math.cos(E)) + r["crs"] * s2 + r["crc"] * c2
    inc = r["i0"] + r["idot"] * tk + r["cis"] * s2 + r["cic"] * c2
    xp, yp = rad * math.cos(u), rad * math.sin(u)
    ci, si = math.cos(inc), math.sin(inc)
    toe = r["toe"]
    if sys_id == "C" and (prn <= 5 or prn >= 59):      # BeiDou GEO
        om = r["om0"] + r["omd"] * tk - we * toe
        co, so = math.cos(om), math.sin(om)
        xg, yg, zg = xp * co - yp * ci * so, xp * so + yp * ci * co, yp * si
        sino, coso = math.sin(we * tk), math.cos(we * tk)
        return (xg * coso + yg * sino * COS_5 + zg * sino * SIN_5,
                -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5,
                -yg * SIN_5 + zg * COS_5)
    om = r["om0"] + (r["omd"] - we) * tk - we * toe
    co, so = math.cos(om), math.sin(om)
    return xp * co - yp * ci * so, xp * so + yp * ci * co, yp * si


def _glo_deriv(s, acc):
    x, y, z, vx, vy, vz = s
    r2 = x * x + y * y + z * z
    r3 = r2 * math.sqrt(r2)
    w2 = GLO_OMEGA ** 2
    a = 1.5 * GLO_J2 * GLO_MU * GLO_RE ** 2 / r2 / r3
    b = 5.0 * z * z / r2
    c = -GLO_MU / r3 - a * (1.0 - b)
    return (vx, vy, vz,
            (c + w2) * x + 2.0 * GLO_OMEGA * vy + acc[0],
            (c + w2) * y - 2.0 * GLO_OMEGA * vx + acc[1],
            (c - 2.0 * a) * z + acc[2])


def glonass_position(r, t):
    """ECEF position (m) of a GLONASS satellite at GPS time t (RK4 integration)."""
    s = tuple(r["pos"]) + tuple(r["vel"])
    acc = r["acc"]
    dt_total = t - r["t"]
    steps = max(1, int(math.ceil(abs(dt_total) / 60.0)))
    h = dt_total / steps
    for _ in range(steps):
        k1 = _glo_deriv(s, acc)
        k2 = _glo_deriv(tuple(s[i] + k1[i] * h / 2 for i in range(6)), acc)
        k3 = _glo_deriv(tuple(s[i] + k2[i] * h / 2 for i in range(6)), acc)
        k4 = _glo_deriv(tuple(s[i] + k3[i] * h for i in range(6)), acc)
        s = tuple(s[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6))
    return s[0], s[1], s[2]


def nearest_record(recs, t, max_dt):
    best = None
    for r in recs:
        d = abs(t - r["t"])
        if d <= max_dt and (best is None or d < best[0]):
            best = (d, r)
    return best[1] if best else None


# ---------------------------------------------------------------------------
# DOP
# ---------------------------------------------------------------------------

def invert(m):
    """Invert a small symmetric matrix (Gauss-Jordan). Returns None if singular."""
    n = len(m)
    a = [list(row) + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(m)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-12:
            return None
        a[col], a[piv] = a[piv], a[col]
        p = a[col][col]
        a[col] = [v / p for v in a[col]]
        for r in range(n):
            if r != col and a[r][col]:
                f = a[r][col]
                a[r] = [vr - f * vc for vr, vc in zip(a[r], a[col])]
    return [row[n:] for row in a]


def dops(los, systems):
    """PDOP, HDOP, VDOP from ENU unit vectors, one receiver clock per system."""
    used = sorted({s for s, _ in los})
    if len(los) < 3 + len(used):
        return None
    rows = [[e, n, u] + [1.0 if s == k else 0.0 for k in used] for s, (e, n, u) in los]
    k = len(rows[0])
    hth = [[sum(r[i] * r[j] for r in rows) for j in range(k)] for i in range(k)]
    q = invert(hth)
    if q is None or min(q[0][0], q[1][1], q[2][2]) < 0:
        return None
    return (math.sqrt(q[0][0] + q[1][1] + q[2][2]),
            math.sqrt(q[0][0] + q[1][1]),
            math.sqrt(q[2][2]))


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def site_location(session):
    proj = SITE_PROJECTION.get(session.get("site"))
    if not proj:
        return None
    lat, lon = lcc_inverse(session["mean_m"]["n"], session["mean_m"]["e"], proj)
    return lat, lon, session["mean_m"]["u"]


def build_geometry(session):
    loc = site_location(session)
    if loc is None:
        return None, "no site tag, so its location is not known"
    lat, lon, h = loc
    rx = geodetic_to_ecef(lat, lon, h)
    rot = enu_rotation(lat, lon)

    start_utc = local_to_utc(session["start"], DATA_TIME_ZONE)
    end_utc = local_to_utc(session["end"], DATA_TIME_ZONE)

    # Orbits for every UTC day touched, plus the neighbours (ephemerides from
    # just before midnight are needed at the start of a day).
    records, nav_files = {}, []
    day = (start_utc - timedelta(hours=4)).date()
    last = (end_utc + timedelta(hours=4)).date()
    while day <= last:
        path = nav_file_for(day)
        nav_files.append(path.name)
        for sat, recs in parse_nav(path).items():
            records.setdefault(sat, []).extend(recs)
        day += timedelta(days=1)

    test_systems = set(session.get("test") or "")
    t0 = utc_to_gps_seconds(start_utc)
    total = (end_utc - start_utc).total_seconds()
    mask = math.radians(ELEV_MASK_DEG)

    out_t = []
    counts = {s: [] for s in SYSTEMS}
    dop_p, dop_h, dop_v = [], [], []
    step = 0
    while step * STEP_S <= total + 1e-6:
        dt = step * STEP_S
        t = t0 + dt
        n_by_sys = {s: 0 for s in SYSTEMS}
        los = []
        for sat, recs in records.items():
            sys_id, prn = sat[0], int(sat[1:])
            r = nearest_record(recs, t, MAX_DTOE[sys_id])
            if r is None or not r["healthy"]:
                continue
            x, y, z = glonass_position(r, t) if sys_id == "R" else kepler_position(sys_id, prn, r, t)
            d = (x - rx[0], y - rx[1], z - rx[2])
            rng = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
            enu = [sum(rot[i][j] * d[j] for j in range(3)) / rng for i in range(3)]
            if math.asin(max(-1.0, min(1.0, enu[2]))) < mask:
                continue
            n_by_sys[sys_id] += 1
            if sys_id in test_systems:
                los.append((sys_id, enu))
        dd = dops(los, test_systems)
        out_t.append(int(dt))
        for s in SYSTEMS:
            counts[s].append(n_by_sys[s])
        for arr, k in ((dop_p, 0), (dop_h, 1), (dop_v, 2)):
            arr.append(None if dd is None else int(round(dd[k] * 100)))
        step += 1

    geometry = {
        "id": session["id"],
        "format": FORMAT_VERSION,
        "start": session["start"],
        "time_zone": DATA_TIME_ZONE,
        "step_s": STEP_S,
        "mask_deg": ELEV_MASK_DEG,
        "location": {"lat": round(lat, 5), "lon": round(lon, 5), "h_m": round(h, 1)},
        "test": session.get("test"),
        "systems": SYSTEMS,
        "nav_files": nav_files,
        "t": out_t,
        "count": counts,
        "dop_scale": 100,
        "pdop": dop_p,
        "hdop": dop_h,
        "vdop": dop_v,
    }
    return geometry, None


def main(argv):
    rebuild_all = "--all" in argv
    only = [a for a in argv if not a.startswith("--")]
    if not INDEX.exists():
        print("No sessions yet.")
        return 0
    sessions = json.loads(INDEX.read_text()).get("sessions", [])
    GEOMETRY.mkdir(parents=True, exist_ok=True)
    failures = 0
    for s in sessions:
        if only and s["id"] not in only:
            continue
        path = GEOMETRY / f"{s['id']}.json"
        if path.exists() and not rebuild_all and not only:
            try:
                old = json.loads(path.read_text())
                if old.get("format") == FORMAT_VERSION and old.get("test") == s.get("test"):
                    continue
            except ValueError:
                pass
        try:
            geometry, reason = build_geometry(s)
        except Exception as err:                   # noqa: BLE001 - report and keep going
            print(f"ERROR computing satellite geometry for {s['id']}: {err}", file=sys.stderr)
            failures += 1
            continue
        if geometry is None:
            print(f"Skipping {s['id']}: {reason}")
            continue
        path.write_text(json.dumps(geometry, separators=(",", ":")))
        p = [v for v in geometry["pdop"] if v is not None]
        print(f"Satellite geometry for {s['id']} ({s.get('site')} / {s.get('test')}): "
              f"{len(geometry['t'])} epochs, at {geometry['location']['lat']:.4f}, "
              f"{geometry['location']['lon']:.4f}; mean visible "
              + ", ".join(f"{k} {sum(v) / len(v):.1f}" for k, v in geometry["count"].items())
              + (f"; PDOP {min(p) / 100:.2f}-{max(p) / 100:.2f}" if p else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
