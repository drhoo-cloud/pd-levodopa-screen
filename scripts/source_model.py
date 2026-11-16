#!/usr/bin/env python3
"""Does isolation source shift tyrDC carriage once species is held fixed?

Usage:  python3 scripts/source_model.py [data dir] [out tsv]
        defaults: data/  and  source_model_stats.tsv

The specification is the one used for the two gates in Table S4b: a logistic
model with species fixed effects, fitted to the species that contain both
outcomes. Here the predictor is the isolation-source category rather than the
glycoside hydrolase count, so the question is whether source adds anything
after lineage is accounted for. It does not, which is why the paper reports
no by-source prevalence.

Needs numpy, pandas, statsmodels and scipy; verify_all.py checks the headline
values from the written table and needs none of them.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

SRC = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "source_model_stats.tsv"

df = pd.read_csv(os.path.join(SRC, "census_isolation_source.tsv"), sep="\t")
df["y"] = df["tyrdc"].astype(int)
df = df.rename(columns={"src_cat": "src"})
res = []

# -- 1. pooled, no species term: the confounded view -----------------------
D = pd.get_dummies(df.src, drop_first=True).astype(float)
ll0 = sm.GLM(df.y, np.ones((len(df), 1)), family=sm.families.Binomial()).fit().llf
ll1 = sm.GLM(df.y, sm.add_constant(D), family=sm.families.Binomial()).fit().llf
lr1 = 2 * (ll1 - ll0)
res.append(["Logistic, source only, no species term", "genome", len(df),
            f"LRT chi2 {lr1:.1f}", f"df {D.shape[1]}",
            f"{stats.chi2.sf(lr1, D.shape[1]):.3g}"])
print(f"[1] pooled           LRT {lr1:,.1f} df={D.shape[1]} "
      f"P={stats.chi2.sf(lr1, D.shape[1]):.3g}   n={len(df):,}")

# -- 2. species fixed effects, mixed species only --------------------------
g = df.groupby("species").y.agg(["sum", "count"])
mixed = g[(g["sum"] > 0) & (g["sum"] < g["count"])].index.tolist()
dm = df[df.species.isin(mixed)].copy()
S = pd.get_dummies(dm.species, drop_first=True).astype(float).reset_index(drop=True)
C = pd.get_dummies(dm.src, drop_first=True).astype(float).reset_index(drop=True)
y = dm.y.reset_index(drop=True)
base = sm.GLM(y, sm.add_constant(S), family=sm.families.Binomial()).fit()
full = sm.GLM(y, sm.add_constant(pd.concat([S, C], axis=1)),
              family=sm.families.Binomial()).fit()
lr2 = 2 * (full.llf - base.llf)
p2 = stats.chi2.sf(lr2, C.shape[1])
res.append(["Logistic, source with species fixed effects", "genome", len(dm),
            f"LRT chi2 {lr2:.2f}", f"df {C.shape[1]}", f"{p2:.4f}"])
print(f"[2] species fixed    LRT {lr2:.2f} df={C.shape[1]} P={p2:.4f}   "
      f"n={len(dm):,}, {len(mixed)} mixed species")
for k in C.columns:
    b, se, p = full.params[k], full.bse[k], full.pvalues[k]
    res.append([f"  source term: {k}", "genome", len(dm),
                f"OR {math.exp(b):.2f}",
                f"{math.exp(b-1.96*se):.2f}-{math.exp(b+1.96*se):.2f}",
                f"{p:.3f}"])

# -- 3. per-species tests live in source_by_species.py ---------------------
# They are not computed here. Several species x source cells hold fewer than
# five genomes, where the chi-square approximation is not safe: an earlier
# version of this script used it, dropped the sources holding fewer than 20
# genomes, and called three of the six species significant. The exact test
# in source_by_species.py permutes the source labels within each species and
# calls one. That script writes census_source_by_species.tsv, which is what
# Table S4c reports and what verify_all.py reads.
res.append(["Per-species tests", "species", len(mixed),
            "see census_source_by_species.tsv",
            "exact, permutation (source_by_species.py)", "—"])

pd.DataFrame(res, columns=["Model or test", "Unit", "n", "Estimate",
                           "95% CI / LRT", "P"]).to_csv(
    OUT, sep="\t", index=False, lineterminator="\n")
print("wrote", OUT)
