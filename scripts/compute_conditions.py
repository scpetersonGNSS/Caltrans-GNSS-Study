#!/usr/bin/env python3
"""
Add satellite and ionosphere conditions to each processed session.

For every tagged session in data/index.json this computes, every STEP_S
seconds through the session:

  Satellites in view   number of satellites above ELEVATION_MASK_DEG at the
                       site, counting only the constellations in that test
                       (GPS only for test G, all four for GREC, ...)
  Ionosphere VTEC      vertical total electron content above the site, in
                       TEC units (1 TECU delays the GPS L1 signal by ~16 cm)

These are the same quantities Trimble GNSS Planning charts, but computed from
the precise orbits and ionosphere maps published by the International GNSS
Service (IGS) for the actual dates of the session, so nothing has to be
copied by hand. The values are stored in the session file next to PDOP, so
the website can plot them the same way.

Data sources (public, no account needed), tried in order:
  orbits      multi-GNSS SP3 files (CODE, GFZ, Wuhan, ...) from the BKG IGS
              mirror and the CODE archive at the University of Bern
  ionosphere  IONEX global ionosphere maps (CODE, IGS) from the same places

Final products appear about one to two weeks after the date and rapid ones
about a day after. Until a session's products exist it is marked "pending"
and the website shows the session without these values; the daily run of
the workflow fills them in once the files are published.

Collector times are local clock times (see COLLECTOR_TIME_ZONE); they are
converted to UTC to match the orbit and ionosphere files.

Needs numpy. Run locally:  python scripts/compute_conditions.py
"""

import argparse
import gzip
import json
import math
import re
import sys
import time
import urllib.error
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

# ---------------------------------------------------------------------------

GPS_EPOCH = date(1980, 1, 6)
GPS_MINUS_UTC = 18                         # seconds, since 2017
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563
USER_AGENT = "Caltrans-GNSS-Study site builder (github.com/scpetersonGNSS/Caltrans-GNSS-Study)"

# Where to look, in order. Each entry is a directory URL and a file-name
# pattern; {y} = year, {d} = day of year (3 digits), {w} = GPS week.
ORBIT_DIRS = [
    "https://igs.bkg.bund.de/root_ftp/IGS/products/mgex/{w}/",
    "http://ftp.aiub.unibe.ch/CODE_MGEX/CODE/{y}/",
    "https://igs.bkg.bund.de/root_ftp/IGS/products/{w}/",
]
ORBIT_RE = r"(?P<ac>[A-Z0-9]{{3}})0(?P<camp>MGX|OPS)(?P<kind>FIN|RAP)_{y}{d}0000_01D_\d\dM_ORB\.SP3\.gz"
ORBIT_GUESSES = [
    "http://ftp.aiub.unibe.ch/CODE_MGEX/CODE/{y}/COD0MGXFIN_{y}{d}0000_01D_05M_ORB.SP3.gz",
]
ION_DIRS = [
    "https://igs.bkg.bund.de/root_ftp/IGS/products/ionex/{y}/{d}/",
    "http://ftp.aiub.unibe.ch/CODE/{y}/",
    "http://ftp.aiub.unibe.ch/CODE/",
]
ION_RE = r"(?P<ac>[A-Z0-9]{{3}})0(?P<camp>OPS|MGX)(?P<kind>FIN|RAP)_{y}{d}0000_01D_\d\d[HM]_GIM\.INX\.gz"
ION_GUESSES = [
    "http://ftp.aiub.unibe.ch/CODE/{y}/COD0OPSFIN_{y}{d}0000_01D_01H_GIM.INX.gz",
    "http://ftp.aiub.unibe.ch/CODE/COD0OPSRAP_{y}{d}0000_01D_01H_GIM.INX.gz",
]
# Preferred analysis centers, best first. Final products beat rapid ones.
ORBIT_CENTERS = ["COD", "GFZ", "WUM", "GRG", "IGS", "ESA", "JAX", "IAC", "SHA"]
ION_CENTERS = ["COD", "IGS", "JPL", "ESA", "UPC", "EMR", "WHU", "CAS"]


# ---------------------------------------------------------------------------
# Downloading
# ---------------------------------------------------------------------------

def fetch(url, tries=3):
    """Return the bytes at url, or None if it does not exist."""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read()
        except urllib.error.HTTPError as err:
            if err.code in (403, 404, 410):
                return None
            last = err
        except (urllib.error.URLError, TimeoutError, OSError) as err:
            last = err
        time.sleep(3 * (attempt + 1))
    print(f"  could not reach {url}: {last}", file=sys.stderr)
    return None


def fmt(template, day):
    week = (day - GPS_EPOCH).days // 7
    return template.format(y=f"{day.year:04d}", d=f"{day.timetuple().tm_yday:03d}", w=f"{week:04d}")


def candidates(day, dirs, pattern, guesses, centers):
    """Product URLs for one UTC day, best first."""
    found = {}
    rx = re.compile(fmt(pattern, day))
    for d in dirs:
        listing = fetch(fmt(d, day))
        if not listing:
            continue
        for m in rx.finditer(listing.decode("utf-8", "replace")):
            found.setdefault(m.group(0), (fmt(d, day) + m.group(0), m))
    def rank(item):
        name, (url, m) = item
        ac = m.group("ac")
        return (m.group("kind") != "FIN", centers.index(ac) if ac in centers else len(centers), m.group("camp") != "MGX", name)
    urls = [url for _, (url, _) in sorted(found.items(), key=rank)]
    for g in guesses:
        u = fmt(g, day)
        if all(Path(u).name != Path(x).name for x in urls):
            urls.append(u)
    return urls


def load_product(url):
    """Download (or reuse) a gzipped product; return its text, or None."""
    CACHE.mkdir(exist_ok=True)
    local = CACHE / Path(url).name
    if local.exists():
        raw = local.read_bytes()
    else:
        raw = fetch(url)
        if raw is None:
            return None
        local.write_bytes(raw)
    try:
        return gzip.decompress(raw).decode("ascii", "replace")
    except (OSError, EOFError):
        local.unlink(missing_ok=True)
        return None


# ---------------------------------------------------------------------------
# Orbits (SP3)
# ---------------------------------------------------------------------------

def parse_sp3(text):
    """{'G01': (times[s, UTC posix], xyz[m] Nx3), ...}"""
    sats, t = {}, None
    system_gps = True
    for line in text.splitlines():
        if line.startswith("%c") and len(line) > 13:
            system_gps = line[9:12].strip() in ("GPS", "ccc", "")
        elif line.startswith("*"):
            p = line[1:].split()
            y, mo, d, h, mi = map(int, p[:5])
            sec = float(p[5])
            t = datetime(y, mo, d, h, mi, tzinfo=timezone.utc).timestamp() + sec
            if system_gps:
                t -= GPS_MINUS_UTC
        elif line.startswith("P") and t is not None and len(line) >= 46:
            prn = line[1:4].replace(" ", "0")
            try:
                xyz = (float(line[4:18]), float(line[18:32]), float(line[32:46]))
            except ValueError:
                continue
            if xyz == (0.0, 0.0, 0.0) or any(abs(v) > 99999 for v in xyz):
                continue                    # no position at this epoch
            sats.setdefault(prn, []).append((t, xyz[0] * 1000, xyz[1] * 1000, xyz[2] * 1000))
    out = {}
    for prn, rows in sats.items():
        a = np.array(sorted(rows))
        out[prn] = (a[:, 0], a[:, 1:4])
    return out


def merge_orbits(parts):
    merged = {}
    for orb in parts:
        for prn, (t, xyz) in orb.items():
            if prn in merged:
                t0, x0 = merged[prn]
                t = np.concatenate([t0, t]); xyz = np.vstack([x0, xyz])
            order = np.argsort(t, kind="stable")
            t, xyz = t[order], xyz[order]
            keep = np.concatenate([[True], np.diff(t) > 0.5])
            merged[prn] = (t[keep], xyz[keep])
    return merged


def sat_positions(orbits, tq):
    """Positions at query times (linear interpolation of 5-15 min samples,
    which is far more precise than satellite counting needs). NaN where a
    satellite has no orbit or the samples are more than 20 minutes apart."""
    out = {}
    for prn, (t, xyz) in orbits.items():
        if len(t) < 2:
            continue
        idx = np.searchsorted(t, tq)
        ok = (idx > 0) & (idx < len(t))
        i1 = np.clip(idx, 1, len(t) - 1)
        ok &= (t[i1] - t[i1 - 1]) <= 1200
        exact = np.isin(tq, t)
        p = np.column_stack([np.interp(tq, t, xyz[:, k]) for k in range(3)])
        p[~(ok | exact)] = np.nan
        out[prn] = p
    return out


def station_ecef(lat, lon, h):
    la, lo = math.radians(lat), math.radians(lon)
    e2 = WGS84_F * (2 - WGS84_F)
    N = WGS84_A / math.sqrt(1 - e2 * math.sin(la) ** 2)
    r = np.array([(N + h) * math.cos(la) * math.cos(lo), (N + h) * math.cos(la) * math.sin(lo), (N * (1 - e2) + h) * math.sin(la)])
    east = np.array([-math.sin(lo), math.cos(lo), 0.0])
    north = np.array([-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la)])
    up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
    return r, np.vstack([east, north, up])


def geometry(orbits, tq, site, letters):
    """Satellites in view (count) and PDOP for the given constellations."""
    r0, enu = station_ecef(*site)
    pos = sat_positions({p: v for p, v in orbits.items() if p[0] in letters}, tq)
    sys_list = sorted(set(p[0] for p in pos))
    count = np.zeros(len(tq), dtype=int)
    rows = []                                   # per satellite: (visible mask, unit los ENU, system index)
    for prn, p in pos.items():
        los = p - r0
        rng = np.linalg.norm(los, axis=1)
        u = (los @ enu.T) / rng[:, None]
        el = np.degrees(np.arcsin(np.clip(u[:, 2], -1, 1)))
        vis = np.nan_to_num(el, nan=-90) >= ELEVATION_MASK_DEG
        count += vis
        rows.append((vis, u, sys_list.index(prn[0])))
    pdop = np.full(len(tq), np.nan)
    ns = len(sys_list)
    for k in range(len(tq)):
        G = [list(-u[k]) + [1.0 if s == j else 0.0 for j in range(ns)] for vis, u, s in rows if vis[k]]
        if len(G) < 3 + ns:
            continue
        G = np.array(G)
        used = G[:, 3:].any(axis=0)             # drop clock columns of systems with no satellites
        G = np.hstack([G[:, :3], G[:, 3:][:, used]])
        try:
            Q = np.linalg.inv(G.T @ G)
        except np.linalg.LinAlgError:
            continue
        pdop[k] = math.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])
    return count, pdop, sys_list


# ---------------------------------------------------------------------------
# Ionosphere (IONEX)
# ---------------------------------------------------------------------------

def parse_ionex(text):
    """List of (time UTC posix, lats, lons, grid[lat, lon] in TECU)."""
    head, _, body = text.partition("END OF HEADER")
    exponent = -1
    lat1 = lat2 = dlat = lon1 = lon2 = dlon = None
    for line in head.splitlines():
        label = line[60:].strip()
        if label == "EXPONENT":
            exponent = int(line[:6])
        elif label == "LAT1 / LAT2 / DLAT":
            lat1, lat2, dlat = map(float, line[:60].split()[:3])
        elif label == "LON1 / LON2 / DLON":
            lon1, lon2, dlon = map(float, line[:60].split()[:3])
    lats = np.arange(lat1, lat2 + dlat / 2, dlat) if dlat else None
    lons = np.arange(lon1, lon2 + dlon / 2, dlon) if dlon else None
    nlon = len(lons)
    maps, grid, epoch, scale = [], None, None, 10.0 ** exponent
    lines = body.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        label = line[60:].strip()
        if label == "START OF RMS MAP" or label == "START OF HEIGHT MAP":
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
                row = lines[i]
                vals += [int(row[k:k + 5]) for k in range(0, len(row.rstrip()), 5) if row[k:k + 5].strip()]
            j = int(round((lat - lat1) / dlat))
            v = np.array(vals[:nlon], dtype=float)
            v[v == 9999] = np.nan
            if 0 <= j < len(lats):
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
        j, i = int(math.floor(fy)), int(math.floor(fx))
        j = min(max(j, 0), len(lats) - 2); i = min(max(i, 0), len(lons) - 2)
        wy, wx = fy - j, fx - i
        return ((1 - wy) * ((1 - wx) * g[j, i] + wx * g[j, i + 1]) + wy * ((1 - wx) * g[j + 1, i] + wx * g[j + 1, i + 1]))
    vals = np.array([at_map(m) for m in maps])
    out = np.interp(tq, times, vals, left=np.nan, right=np.nan)
    gap = np.searchsorted(times, tq)
    gi = np.clip(gap, 1, len(times) - 1)
    out[(times[gi] - times[gi - 1]) > 3 * 3600 + 1] = np.nan
    return out


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def utc_days(t0, t1):
    d, last = (datetime.fromtimestamp(t0, timezone.utc) - timedelta(minutes=30)).date(), \
              (datetime.fromtimestamp(t1, timezone.utc) + timedelta(minutes=30)).date()
    days = []
    while d <= last:
        days.append(d); d += timedelta(days=1)
    return days


def get_orbits(days, letters, log):
    parts, used = [], []
    for day in days:
        got = None
        for url in candidates(day, ORBIT_DIRS, ORBIT_RE, ORBIT_GUESSES, ORBIT_CENTERS):
            text = load_product(url)
            if not text:
                continue
            orb = parse_sp3(text)
            have = {p[0] for p in orb}
            if not set(letters) <= have:
                log(f"  {Path(url).name} lacks {''.join(sorted(set(letters) - have))}; trying next")
                continue
            got = (orb, Path(url).name)
            break
        if not got:
            return None, f"no orbit file with {letters} yet for {day}"
        parts.append(got[0]); used.append(got[1])
    return (merge_orbits(parts), used), None


def get_ionosphere(days, log):
    maps, used = [], []
    for day in days:
        got = None
        for url in candidates(day, ION_DIRS, ION_RE, ION_GUESSES, ION_CENTERS):
            text = load_product(url)
            if not text:
                continue
            m = parse_ionex(text)
            if m:
                got = (m, Path(url).name)
                break
        if not got:
            return None, f"no ionosphere map yet for {day}"
        maps += got[0]; used.append(got[1])
    return (maps, used), None


def session_times(data, offset_hours=None):
    """UTC posix time of the session start, and the computation grid (s from start)."""
    start = datetime.fromisoformat(data["start"])
    if offset_hours is None:
        start_utc = start.replace(tzinfo=ZoneInfo(COLLECTOR_TIME_ZONE)).timestamp()
    else:
        start_utc = start.replace(tzinfo=timezone(timedelta(hours=offset_hours))).timestamp()
    end = data["t"][-1]
    grid = np.arange(0, end + STEP_S, STEP_S, dtype=float)
    grid = grid[grid <= end + STEP_S / 2]
    return start_utc, grid


def as_extra(values, dec, unit, note):
    scale = 10 ** dec
    return {"unit": unit, "dec": dec, "scale": scale, "step_s": STEP_S, "computed": True, "note": note,
            "v": [None if not np.isfinite(v) else int(round(v * scale)) for v in values]}


def compute(summary, data, log):
    site = SITE_LOCATIONS[summary["site"]]
    letters = summary["test"]
    start_utc, grid = session_times(data)
    tq = start_utc + grid
    days = utc_days(tq[0], tq[-1])
    orb, why = get_orbits(days, letters, log)
    if not orb:
        return None, why
    ion, why = get_ionosphere(days, log)
    if not ion:
        return None, why
    (orbits, orb_files), (maps, ion_files) = orb, ion
    count, _, _ = geometry(orbits, tq, site, letters)
    vtec = vtec_at(maps, tq, site[0], site[1])
    cons = " + ".join(summary.get("constellations", [])) or letters
    extras = {
        SAT_NAME: as_extra(count.astype(float), 0, "",
                           f"Computed: {cons} satellites above {ELEVATION_MASK_DEG:g}° from IGS precise orbits ({', '.join(orb_files)})."),
        ION_NAME: as_extra(vtec, 1, "TECU",
                           f"Computed: vertical TEC above the site from IGS ionosphere maps ({', '.join(ion_files)}). 1 TECU ≈ 16 cm of L1 delay."),
    }
    info = {"status": "done", "version": CONDITIONS_VERSION, "mask_deg": ELEVATION_MASK_DEG,
            "step_s": STEP_S, "time_zone": COLLECTOR_TIME_ZONE, "orbits": orb_files, "ionosphere": ion_files,
            "computed": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    return (extras, info, orbits, tq), None


def recorded_on_grid(data, grid):
    """Recorded PDOP sampled at the grid times (nearest epoch within 5 s)."""
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


def time_zone_check(summary, data, orbits):
    """Correlation of computed vs recorded PDOP for a range of UTC offsets."""
    site = SITE_LOCATIONS[summary["site"]]
    rows = []
    for off in range(-10, 3):
        start_utc, grid = session_times(data, off)
        _, pdop, _ = geometry(orbits, start_utc + grid, site, summary["test"])
        rec = recorded_on_grid(data, grid)
        if rec is None:
            return None
        ok = np.isfinite(pdop) & np.isfinite(rec)
        r = float(np.corrcoef(pdop[ok], rec[ok])[0, 1]) if ok.sum() > 10 else float("nan")
        rows.append((off, r, int(ok.sum())))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", help="comma-separated session ids to compute")
    ap.add_argument("--force", action="store_true", help="recompute even if already done")
    ap.add_argument("--report", help="write a check report (Markdown) to this file")
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
        print(f"{s['id']} ({s['site']} / {s['test']}):")
        result, why = compute(s, data, print)
        if not result:
            print(f"  pending: {why}")
            new = {"status": "pending", "reason": why, "version": CONDITIONS_VERSION,
                   "checked": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
            if cond.get("status") != "done":
                changed |= cond != new
                s["conditions"] = new
            report += [f"## {s['id']} {s['site']} {s['test']}", f"- pending: {why}", ""]
            continue
        extras, info, orbits, tq = result
        sats = [v for v in extras[SAT_NAME]["v"] if v is not None]
        ions = [v / 10 for v in extras[ION_NAME]["v"] if v is not None]
        print(f"  satellites {min(sats)}–{max(sats)}, VTEC {min(ions):.1f}–{max(ions):.1f} TECU "
              f"from {', '.join(info['orbits'])} and {', '.join(info['ionosphere'])}")
        if args.report:
            report += [f"## {s['id']} {s['site']} {s['test']}",
                       f"- orbits: {', '.join(info['orbits'])}", f"- ionosphere: {', '.join(info['ionosphere'])}",
                       f"- satellites: min {min(sats)}, max {max(sats)}, mean {sum(sats)/len(sats):.1f}",
                       f"- VTEC TECU: min {min(ions):.1f}, max {max(ions):.1f}",
                       "- first hour (local, every 10 min): " + ", ".join(
                           f"{extras[SAT_NAME]['v'][k]}/{extras[ION_NAME]['v'][k]/10 if extras[ION_NAME]['v'][k] is not None else '-'}"
                           for k in range(0, min(121, len(extras[SAT_NAME]['v'])), 20))]
            rows = time_zone_check(s, data, orbits)
            if rows:
                best = max(rows, key=lambda r: -1 if math.isnan(r[1]) else r[1])
                report.append(f"- PDOP correlation by UTC offset (best {best[0]:+d} h): " +
                              ", ".join(f"{o:+d}h {r:.2f}" for o, r, _ in rows))
            report.append("")
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
