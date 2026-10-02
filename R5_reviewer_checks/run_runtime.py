#!/usr/bin/env python3
"""Runtime and memory of one MDFS run per data type (Reviewer 1, comment 9).
LOCO fold-0 training set; run_mdfs.R with seed 42 under /usr/bin/time -l (macOS).
Output: runtime.tsv (samples, features, candidate pairs n(n-1)/2, wall-clock seconds, peak RSS).
"""
import os
import platform
import subprocess
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import training_set, run_mdfs, unordered, HERE  # noqa: E402

rows = []
for dt in ("taxa", "functions", "combined"):
    X, y = training_set(dt)
    pairs, ig1d, out, wall, rss = run_mdfs(X, y, 42)
    rows.append({"data_type": dt, "n_samples": len(y), "n_features": X.shape[1], "candidate_pairs": X.shape[1] * (X.shape[1] - 1) // 2,
                 "pairs_passing_filters": len(unordered(pairs)), "wall_clock_s": round(wall, 1),
                 "peak_RSS_GB": round(rss / 1e9, 2) if rss else None})
    print(rows[-1], flush=True)
cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True).stdout.strip()
ncpu = subprocess.run(["sysctl", "-n", "hw.ncpu"], capture_output=True, text=True).stdout.strip()
mem = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout) / 1e9
rv = subprocess.run(["Rscript", "-e", 'cat(R.version.string, "MDFS", as.character(packageVersion("MDFS")))'], capture_output=True, text=True).stdout
out = pd.DataFrame(rows)
out["machine"] = f"{cpu}, {ncpu} cores, {mem:.0f} GB RAM, {platform.platform()}; {rv}"
out.to_csv(os.path.join(HERE, "runtime.tsv"), sep="\t", index=False)
print(out.to_string(index=False))
