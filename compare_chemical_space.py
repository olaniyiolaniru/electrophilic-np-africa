"""Reproduce the property-matched global chemical-space comparison: the difference
between motif-bearing and matched motif-negative structures in the fraction of
records below each nearest-neighbor Tanimoto threshold to the frozen COCONUT reference.

Inputs (data/):
    global_matched_pairs.tsv        3,416 matched positive/control pairs
    global_nearest_neighbors.tsv    per-record maximum Tanimoto at 2048 and 4096 bits
Output (out/):
    global_comparison.tsv           per (bits, threshold) gap with paired 95% interval

Run:  python compare_chemical_space.py
Requires: numpy, pandas
"""
from pathlib import Path
import numpy as np, pandas as pd

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)
THRESHOLDS = [0.5, 0.6, 0.7, 0.8]; N_BOOT = 2000; SEED = 42


def main():
    pairs = pd.read_csv(DATA / "global_matched_pairs.tsv", sep="\t")
    nn = pd.read_csv(DATA / "global_nearest_neighbors.tsv", sep="\t")
    rows = []
    for bits, gb in nn.groupby("bits"):
        sim = gb.set_index("molecule_id").tanimoto
        pos = pairs.positive_id.map(sim).values; ctl = pairs.control_id.map(sim).values
        rng = np.random.default_rng(SEED); n = len(pos); idx = rng.integers(0, n, size=(N_BOOT, n))
        for thr in THRESHOLDS:
            pb = pos < thr; cb = ctl < thr
            gap = 100 * (pb.mean() - cb.mean())
            boot = 100 * (pb[idx].mean(axis=1) - cb[idx].mean(axis=1))
            lo, hi = np.percentile(boot, [2.5, 97.5])
            rows.append(dict(bits=int(bits), threshold=thr, positive_percent=100 * pb.mean(),
                             control_percent=100 * cb.mean(), gap_pp=gap, ci_low=lo, ci_high=hi))
    out = pd.DataFrame(rows); out.to_csv(OUT / "global_comparison.tsv", sep="\t", index=False)
    for r in out[out.bits == 2048].itertuples():
        print(f"2048 bits, T<{r.threshold}: gap {r.gap_pp:+.2f} pp  [{r.ci_low:.2f}, {r.ci_high:.2f}]")


if __name__ == "__main__":
    main()
