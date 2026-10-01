#!/usr/bin/env python3
"""
Process GNSS static occupation files dropped into uploads/.

Each uploaded file becomes its own session:
  data/sessions/<session-id>.json   plot data (integer mm offsets from a base)
  data/raw/<session-id>_<name>      the original file, for download
  data/index.json                   list of sessions with summary statistics

Every upload must be tagged with its test site and constellation test through
its file name: the site code, then the test code, e.g.

  SJER_G.asc       Site 1 (San Joaquin Experimental Range), GPS only
  FRES_GREC.asc    Site 2 (Fresno State), GPS + GLONASS + Galileo + BeiDou
  ELK_GRC.asc      Site 3 (Elkhorn), GPS + GLONASS + BeiDou

Anything may come before or after the two codes (a date, a note), as long as
the test code directly follows the site code. Letters within a test code may
be in any order (GERC is read as GREC). A file without valid codes is not
processed: the run fails and the file stays in uploads/ so it can be renamed.

Uses only the Python standard library. Column names are auto-detected
(Northing / Easting / Elevation / time), so exports with slightly different
headers still work. UTF-16 exports from the data collector are read directly.

Run locally:   python scripts/process_uploads.py
"""

import csv
import hashlib
import io
import json
import math
import re
import shutil
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UPLOADS = ROOT / "uploads"
SESSIONS = ROOT / "data" / "sessions"
RAW = ROOT / "data" / "raw"
INDEX = ROOT / "data" / "index.json"

ACCEPTED = {".asc", ".csv", ".txt"}

# ---------------------------------------------------------------------------
# Study design. Edit here to add or rename sites and tests.
# ---------------------------------------------------------------------------

# Site code used in file names -> site number and full name.
SITES = {
    "SJER": {"number": 1, "name": "San Joaquin Experimental Range"},
    "FRES": {"number": 2, "name": "Fresno State"},
    "ELK":  {"number": 3, "name": "Elkhorn"},
}

# Constellation letters.
CONSTELLATIONS = {"G": "GPS", "R": "GLONASS", "E": "Galileo", "C": "BeiDou"}

# The five tests run at every site, in order (test 1 first).
TESTS = ["G", "GREC", "GEC", "GRE", "GRC"]

# Sessions processed before tagging existed: session id -> (site, test).
LEGACY_TAGS = {
    "20260926-082350": ("SJER", "G"),
}

# ---------------------------------------------------------------------------

# Timestamp formats seen in collector exports; the first one that parses wins.
TIME_FORMATS = [
    "%m/%d/%Y %I:%M:%S %p",   # 8/31/2026 7:22:47 AM
    "%m/%d/%Y %H:%M:%S",      # 8/31/2026 19:22:47
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%dT%H:%M:%S",
    "%m/%d/%Y %I:%M %p",
    "%m/%d/%Y %H:%M",
]

# A gap longer than this many epoch intervals breaks the plotted line.
GAP_FACTOR = 3

# Bump when the session file layout changes; older sessions are rebuilt
# automatically from their original file in data/raw/.
FORMAT_VERSION = 2

# Summary fields that describe the upload rather than its contents; they are
# kept when a session is rebuilt.
TAG_FIELDS = ("site", "test", "test_number", "constellations")
KEEP_FIELDS = ("id", "source_file", "uploaded", "sha256", "raw_file") + TAG_FIELDS


# ---------------------------------------------------------------------------
# Site / test tags
# ---------------------------------------------------------------------------

def canonical_test(code):
    """Return the test code matching these constellation letters, or None."""
    letters = set(code)
    if len(letters) != len(code) or not letters <= set(CONSTELLATIONS):
        return None
    for test in TESTS:
        if set(test) == letters:
            return test
    return None


def parse_tags(filename):
    """Read (site, test) from a file name such as 'SJER_GREC.asc'."""
    tokens = [t for t in re.split(r"[^A-Za-z0-9]+", Path(filename).stem.upper()) if t]
    for i, token in enumerate(tokens):
        if token in SITES:
            test = canonical_test(tokens[i + 1]) if i + 1 < len(tokens) else None
            if test:
                return token, test
            raise ValueError(
                f"site code {token} found, but it is not followed by a valid test code. "
                f"Tests are: {', '.join(TESTS)}. Example: {token}_GREC.asc")
    raise ValueError(
        "file name must include a site code followed by a test code, "
        f"e.g. SJER_G.asc or FRES_GREC.asc. Sites are: {', '.join(SITES)}. "
        f"Tests are: {', '.join(TESTS)}.")


def tag_fields(site, test):
    return {
        "site": site,
        "test": test,
        "test_number": TESTS.index(test) + 1,
        "constellations": [CONSTELLATIONS[c] for c in test],
    }


def study_metadata():
    """Site and test definitions, written to index.json for the website."""
    return {
        "sites": [{"code": code, **info} for code, info in
                  sorted(SITES.items(), key=lambda kv: kv[1]["number"])],
        "tests": [{"code": t, "number": i + 1,
                   "constellations": [CONSTELLATIONS[c] for c in t]}
                  for i, t in enumerate(TESTS)],
        "constellations": CONSTELLATIONS,
    }


def apply_legacy_tags(sessions):
    changed = False
    for s in sessions:
        if "site" not in s and s["id"] in LEGACY_TAGS:
            s.update(tag_fields(*LEGACY_TAGS[s["id"]]))
            print(f"Tagged existing session {s['id']} as {s['site']} / {s['test']}")
            changed = True
    return changed


# ---------------------------------------------------------------------------
# Reading and summarizing observation files
# ---------------------------------------------------------------------------

def find_column(headers, *keywords):
    """Return the index of the first header containing any keyword (case-insensitive)."""
    lowered = [h.strip().lower() for h in headers]
    for kw in keywords:
        for i, h in enumerate(lowered):
            if kw in h:
                return i
    return None


def parse_time(text):
    text = text.strip()
    for fmt in TIME_FORMATS:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    raise ValueError(f"Unrecognized timestamp: {text!r}")


def read_text(path):
    """Read a text export in whatever encoding the collector used."""
    raw = path.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):       # UTF-16 with byte-order mark
        return raw.decode("utf-16")
    if raw.startswith(b"\xef\xbb\xbf"):                  # UTF-8 with byte-order mark
        return raw[3:].decode("utf-8")
    if len(raw) > 1 and raw[1:2] == b"\x00":             # UTF-16 LE without a mark
        return raw.decode("utf-16-le")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("cp1252", errors="replace")    # older Windows exports


def read_observations(path):
    text = read_text(path)
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    rows = list(csv.reader(io.StringIO(text, newline=""), dialect))

    if not rows:
        raise ValueError("File is empty")
    headers = rows[0]
    col_n = find_column(headers, "northing", "north")
    col_e = find_column(headers, "easting", "east")
    col_u = find_column(headers, "elevation", "height", "elev", "ortho")
    col_t = find_column(headers, "start time", "time", "date")
    col_id = find_column(headers, "point id", "point", "id")
    missing = [name for name, c in
               [("Northing", col_n), ("Easting", col_e), ("Elevation", col_u), ("Time", col_t)]
               if c is None]
    if missing:
        raise ValueError(f"Could not find column(s): {', '.join(missing)}. Headers were: {headers}")

    obs = {}
    skipped = 0
    for row in rows[1:]:
        if not row or all(not c.strip() for c in row):
            continue
        try:
            t = parse_time(row[col_t])
            n = float(row[col_n]); e = float(row[col_e]); u = float(row[col_u])
        except (ValueError, IndexError):
            skipped += 1
            continue
        pid = row[col_id].strip() if col_id is not None and col_id < len(row) else ""
        obs[t] = (n, e, u, pid)   # duplicate timestamps: last one wins
    if not obs:
        raise ValueError("No valid observation rows found")
    times = sorted(obs)
    return times, [obs[t] for t in times], skipped


def std(values):
    return statistics.pstdev(values) if len(values) > 1 else 0.0


def build_session(path):
    times, coords, skipped = read_observations(path)
    N = [c[0] for c in coords]; E = [c[1] for c in coords]; U = [c[2] for c in coords]
    mean_n, mean_e, mean_u = statistics.fmean(N), statistics.fmean(E), statistics.fmean(U)

    # Round-number bases; values are stored as integer millimetres from these,
    # which keeps files small and avoids floating-point transcription errors.
    base = {"n": float(math.floor(mean_n)), "e": float(math.floor(mean_e)), "u": float(math.floor(mean_u))}

    start = times[0]
    t_off = [int(round((t - start).total_seconds())) for t in times]
    diffs = [b - a for a, b in zip(t_off, t_off[1:]) if b > a]
    interval = statistics.median(diffs) if diffs else 0

    sn, se, su = std(N), std(E), std(U)
    drms = math.sqrt(sn ** 2 + se ** 2)
    gaps = sum(1 for d in diffs if interval and d > GAP_FACTOR * interval)

    session_id = start.strftime("%Y%m%d-%H%M%S")
    data = {
        "id": session_id,
        "start": start.strftime("%Y-%m-%dT%H:%M:%S"),
        "interval_s": interval,
        "gap_factor": GAP_FACTOR,
        "base_m": base,
        "t": t_off,
        "n_mm": [int(round((v - base["n"]) * 1000)) for v in N],
        "e_mm": [int(round((v - base["e"]) * 1000)) for v in E],
        "u_mm": [int(round((v - base["u"]) * 1000)) for v in U],
        "ids": [c[3] for c in coords],
        "format": FORMAT_VERSION,
    }
    summary = {
        "id": session_id,
        "source_file": path.name,
        "start": start.strftime("%Y-%m-%dT%H:%M:%S"),
        "end": times[-1].strftime("%Y-%m-%dT%H:%M:%S"),
        "duration_h": round((times[-1] - start).total_seconds() / 3600, 2),
        "epochs": len(times),
        "interval_s": interval,
        "gaps": gaps,
        "skipped_rows": skipped,
        "mean_m": {"n": round(mean_n, 3), "e": round(mean_e, 3), "u": round(mean_u, 3)},
        "std_mm": {"n": round(sn * 1000, 1), "e": round(se * 1000, 1), "u": round(su * 1000, 1)},
        "drms_mm": round(drms * 1000, 1),
        "uploaded": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    return data, summary


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refresh_outdated(sessions):
    """Rebuild session files written by an older version of this script."""
    changed = False
    for s in sessions:
        path = SESSIONS / f"{s['id']}.json"
        raw = ROOT / s.get("raw_file", "")
        try:
            current = json.loads(path.read_text()).get("format", 1)
        except (OSError, ValueError):
            current = 0
        if current >= FORMAT_VERSION or not raw.is_file():
            continue
        try:
            data, summary = build_session(raw)
        except ValueError as err:
            print(f"Could not refresh {s['id']}: {err}", file=sys.stderr)
            continue
        data["id"] = s["id"]
        path.write_text(json.dumps(data, separators=(",", ":")))
        keep = {k: s[k] for k in KEEP_FIELDS if k in s}
        s.clear(); s.update(summary); s.update(keep)
        print(f"Refreshed session {s['id']} to format {FORMAT_VERSION}")
        changed = True
    return changed


# ---------------------------------------------------------------------------

def main():
    index = json.loads(INDEX.read_text()) if INDEX.exists() else {"sessions": []}
    sessions = index.get("sessions", [])
    changed = refresh_outdated(sessions)
    changed = apply_legacy_tags(sessions) or changed

    meta = study_metadata()
    if index.get("study") != meta:
        index["study"] = meta
        changed = True

    known_hashes = {s.get("sha256") for s in sessions}
    used_ids = {s["id"] for s in sessions}

    UPLOADS.mkdir(exist_ok=True)
    uploads = sorted(p for p in UPLOADS.iterdir()
                     if p.is_file() and p.suffix.lower() in ACCEPTED)
    if not uploads and not changed:
        print("No new uploads.")
        return 0

    failures = 0
    for path in uploads:
        digest = file_hash(path)
        if digest in known_hashes:
            print(f"Skipping {path.name}: identical file already processed.")
            path.unlink()
            continue
        try:
            site, test = parse_tags(path.name)
            data, summary = build_session(path)
        except ValueError as err:
            print(f"ERROR in {path.name}: {err}", file=sys.stderr)
            failures += 1
            continue   # leave the file in uploads/ so it can be fixed

        # Two uploads starting in the same second get a suffix.
        sid, k = summary["id"], 2
        while sid in used_ids:
            sid = f"{summary['id']}-{k}"; k += 1
        data["id"] = summary["id"] = sid
        summary["sha256"] = digest
        summary.update(tag_fields(site, test))

        raw_name = f"{sid}_{path.name.replace(' ', '_')}"
        summary["raw_file"] = f"data/raw/{raw_name}"
        (SESSIONS / f"{sid}.json").write_text(json.dumps(data, separators=(",", ":")))
        shutil.move(str(path), RAW / raw_name)

        sessions.append(summary)
        used_ids.add(sid); known_hashes.add(digest)
        print(f"Added session {sid} from {path.name} ({site} / {test}): "
              f"{summary['epochs']} epochs, DRMS {summary['drms_mm']} mm")

    sessions.sort(key=lambda s: s["start"], reverse=True)
    index["sessions"] = sessions
    index["updated"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    INDEX.write_text(json.dumps(index, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
