"""Temporary check: predicted PDOP vs the PDOP logged by the receiver."""
import json, math, statistics, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
out = []
for g in sorted((ROOT / "data" / "geometry").glob("*.json")):
    geo = json.loads(g.read_text())
    ses = json.loads((ROOT / "data" / "sessions" / g.name).read_text())
    ex = ses.get("extra", {}).get("PDOP")
    if not ex:
        continue
    logged = {}
    for t, v in zip(ses["t"], ex["v"]):
        if v is not None:
            logged.setdefault(t // 60, []).append(v / ex["scale"])
    logged = {k: statistics.fmean(v) for k, v in logged.items()}
    pred = {t // 60: p / 100 for t, p in zip(geo["t"], geo["pdop"]) if p is not None}
    line = [g.stem, geo["test"]]
    best = None
    for off in range(-600, 601, 1):          # minutes
        pairs = [(pred[k + off], v) for k, v in logged.items() if k + off in pred]
        if len(pairs) < 300:
            continue
        a, b = zip(*pairs)
        try:
            r = statistics.correlation(a, b)
        except statistics.StatisticsError:
            continue
        if best is None or r > best[1]:
            best = (off, r, statistics.fmean(abs(x - y) for x, y in pairs), statistics.fmean(a), statistics.fmean(b))
    pairs0 = [(pred[k], v) for k, v in logged.items() if k in pred]
    a, b = zip(*pairs0)
    line.append(f"offset0: r={statistics.correlation(a,b):.3f} mean pred {statistics.fmean(a):.2f} logged {statistics.fmean(b):.2f}")
    line.append(f"best offset {best[0]} min r={best[1]:.3f} mad={best[2]:.2f}")
    line.append("counts " + " ".join(f"{k}:{min(v)}-{max(v)}" for k, v in geo["count"].items()))
    out.append(" | ".join(line))
Path(ROOT / "validation.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))
