#!/usr/bin/env python3
"""
Process GNSS static occupation files dropped into uploads/.

Each uploaded file becomes its own session:
  data/sessions/<session-id>.json   plot data (integer mm offsets from a base)
  data/raw/<session-id>_<name>      the original file, for download
  data/index.json                   list of sessions with summary statistics

Uses only the Python standard library. Column names are auto-detected
(Northing / Easting / Elevation / time), so exports with slightly different
headers still work.

Run locally:   python scripts/process_uploads.py
"""

import csv
import hashlib
import io
import json
import math
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
        obs[t] = (n, e, u)   # duplicate timestamps: last one wins
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


def main():
    index = json.loads(INDEX.read_text()) if INDEX.exists() else {"sessions": []}
    sessions = index.get("sessions", [])
    known_hashes = {s.get("sha256") for s in sessions}
    used_ids = {s["id"] for s in sessions}

    uploads = sorted(p for p in UPLOADS.iterdir()
                     if p.is_file() and p.suffix.lower() in ACCEPTED)
    if not uploads:
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

        raw_name = f"{sid}_{path.name.replace(' ', '_')}"
        summary["raw_file"] = f"data/raw/{raw_name}"
        (SESSIONS / f"{sid}.json").write_text(json.dumps(data, separators=(",", ":")))
        shutil.move(str(path), RAW / raw_name)

        sessions.append(summary)
        used_ids.add(sid); known_hashes.add(digest)
        print(f"Added session {sid} from {path.name}: {summary['epochs']} epochs, "
              f"DRMS {summary['drms_mm']} mm")

    sessions.sort(key=lambda s: s["start"], reverse=True)
    index["sessions"] = sessions
    index["updated"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    INDEX.write_text(json.dumps(index, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
