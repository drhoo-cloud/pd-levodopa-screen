#!/usr/bin/env python3
"""Carriage against isolation source, within each species, on the full census.

Usage:  python3 scripts/source_by_species.py [data dir] [out dir]
        defaults: data/  and  the current directory

Writes census_source_by_species.tsv (one row per species) and
census_source_by_species_long.tsv (one row per species x source).

Several species x source cells are small, so the chi-square approximation is
not safe: the P value is exact, built by permuting the source labels within
the species (50,000 permutations, seed 20261112), and corrected by Holm
across the species that contain both outcomes.  Columns restricted to
sources with at least 10 genomes are given alongside the unrestricted ones,
because a category holding one genome gives a 0% or 100% extreme that is not
a spread.
"""
import csv, collections, json, os, sys

import numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else "data"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "."
B, SEED, MIN_N = 50000, 20261112, 10
rng = np.random.default_rng(SEED)

rows = list(csv.DictReader(open(os.path.join(SRC, "census_isolation_source.tsv")),
                           delimiter="\t"))
per = collections.defaultdict(lambda: ([], []))
for r in rows:
    per[r["species"]][0].append(r["src_cat"])
    per[r["species"]][1].append(int(r["tyrdc"]))


def chi2(counts, n_cat, k_tot, n):
    ck = np.stack([counts, n_cat - counts], axis=2)
    cs = np.array([k_tot, n - k_tot], float)
    exp = n_cat[None, :, None] * cs[None, None, :] / n
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(exp > 0, (ck - exp) ** 2 / exp, 0.0).sum(axis=(1, 2))


stat = {}
for sp, (src, y) in per.items():
    y = np.array(y); n = len(y); k = int(y.sum())
    if not (0 < k < n):
        continue
    cats = sorted(set(src)); idx = np.array([cats.index(s) for s in src])
    n_cat = np.bincount(idx, minlength=len(cats)).astype(float)
    obs = np.bincount(idx[y == 1], minlength=len(cats)).astype(float)
    s0 = chi2(obs[None, :], n_cat, k, n)[0]
    hits = 0
    for start in range(0, B, 2000):
        b = min(2000, B - start)
        yp = y[np.argsort(rng.random((b, n)), axis=1)]
        cnt = np.stack([np.bincount(idx[yp[i] == 1], minlength=len(cats))
                        for i in range(b)]).astype(float)
        hits += int((chi2(cnt, n_cat, k, n) >= s0 - 1e-9).sum())
    stat[sp] = dict(chi2=s0, k_sources=len(cats), p=(hits + 1) / (B + 1), n=n)

order = sorted(stat, key=lambda s: stat[s]["p"]); m = len(order); run = 0.0
for i, s in enumerate(order):
    run = max(run, min(1.0, (m - i) * stat[s]["p"])); stat[s]["p_holm"] = run

cell = collections.defaultdict(lambda: [0, 0])
for r in rows:
    c = cell[(r["species"], r["src_cat"])]
    c[0] += int(r["tyrdc"]); c[1] += 1
CATS = sorted({r["src_cat"] for r in rows})

head = ["species", "n_genomes", "n_carriers", "prevalence_pct",
        "n_sources_present", "min_pct", "max_pct", "spread_points",
        f"n_sources_min{MIN_N}", f"min_pct_min{MIN_N}", f"max_pct_min{MIN_N}",
        f"spread_points_min{MIN_N}", "chi2", "P_exact", "P_holm", "tested"]
out = [head]
for sp in sorted(per, key=lambda s: -len(per[s][1])):
    y = per[sp][1]; n = len(y); k = sum(y)
    pres = [(c, cell[(sp, c)]) for c in CATS if cell[(sp, c)][1] > 0]
    pct = [100 * a / b for _, (a, b) in pres]
    big = [100 * a / b for _, (a, b) in pres if b >= MIN_N]
    row = [sp, n, k, "%.2f" % (100 * k / n), len(pres),
           "%.2f" % min(pct), "%.2f" % max(pct), "%.2f" % (max(pct) - min(pct)),
           len(big), "%.2f" % min(big), "%.2f" % max(big),
           "%.2f" % (max(big) - min(big))]
    s = stat.get(sp)
    row += (["%.2f" % s["chi2"], "%.5f" % s["p"], "%.4f" % s["p_holm"], "yes"]
            if s else ["", "", "", "no"])
    out.append(row)
with open(os.path.join(OUTDIR, "census_source_by_species.tsv"), "w", newline="") as f:
    csv.writer(f, delimiter="\t", lineterminator="\n").writerows(out)

out2 = [["species", "src_cat", "genomes", "carriers", "prevalence_pct"]]
for sp in sorted(per, key=lambda s: -len(per[s][1])):
    for c in CATS:
        a, b = cell[(sp, c)]
        if b:
            out2.append([sp, c, b, a, "%.2f" % (100 * a / b)])
with open(os.path.join(OUTDIR, "census_source_by_species_long.tsv"), "w", newline="") as f:
    csv.writer(f, delimiter="\t", lineterminator="\n").writerows(out2)

json.dump(stat, open(os.path.join(OUTDIR, "source_by_species_stats.json"), "w"), indent=1)
print("species %d | tested %d | significant after Holm %d"
      % (len(per), m, sum(1 for s in stat if stat[s]["p_holm"] < 0.05)))
for s in order:
    d = stat[s]
    print(f"  {s:32s} n {d['n']:5d}  chi2 {d['chi2']:7.2f}  P {d['p']:.5f}  Holm {d['p_holm']:.4f}")
