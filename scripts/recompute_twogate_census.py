#!/usr/bin/env python3
"""Recompute Table S4b on the full census.

Input   gate1_gh_families_per_genome.tsv (19,765 genomes, gh_families + tyrdc)
        rpoB_species_tree.nwk, 16S_species_tree.nwk
Output  twogate_census_stats.tsv
"""
import csv, os, sys, math, collections
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import minimize_scalar

SRC = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "twogate_census_stats.tsv"

rows = list(csv.DictReader(open(os.path.join(SRC, "gate1_gh_census_per_genome.tsv"),
                                errors="replace"), delimiter="\t"))
df = pd.DataFrame({"species": [r["species"] for r in rows],
                   "gh": [int(r["gh_families"]) for r in rows],
                   "y":  [1 if r["tyrdc"] == "1" else 0 for r in rows]})
print(f"genomes {len(df):,}  carriers {df.y.sum():,}  species {df.species.nunique()}")

res = []

# ── 1. logistic, no species term -------------------------------------------
X = sm.add_constant(df[["gh"]].astype(float))
m1 = sm.GLM(df.y, X, family=sm.families.Binomial()).fit()
b, se = m1.params["gh"], m1.bse["gh"]
or1, lo1, hi1 = math.exp(b), math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)
z = b / se
logp = stats.norm.logsf(abs(z)) + math.log(2)
res.append(["Logistic, no species term", "genome", len(df),
            f"OR {or1:.3f} per GH family", f"{lo1:.3f}-{hi1:.3f}",
            f"~10^{logp/math.log(10):.0f}"])
print(f"[1] logistic, no species term: OR {or1:.3f} ({lo1:.3f}-{hi1:.3f}), "
      f"P ~ 10^{logp/math.log(10):.0f}   n={len(df):,}")

# -- 2. logistic with species fixed effects (mixed species only) -----------
g = df.groupby("species").y.agg(["sum", "count"])
mixed = g[(g["sum"] > 0) & (g["sum"] < g["count"])].index.tolist()
dm = df[df.species.isin(mixed)].copy()
print(f"[2] mixed species {len(mixed)}, genomes {len(dm):,}")
D = pd.get_dummies(dm.species, drop_first=True).astype(float)
X2 = sm.add_constant(pd.concat([dm[["gh"]].astype(float).reset_index(drop=True),
                                D.reset_index(drop=True)], axis=1))
m2 = sm.GLM(dm.y.reset_index(drop=True), X2, family=sm.families.Binomial()).fit()
b2, se2, p2 = m2.params["gh"], m2.bse["gh"], m2.pvalues["gh"]
res.append(["Logistic, species fixed effects", "genome", len(dm),
            f"OR {math.exp(b2):.3f}",
            f"{math.exp(b2-1.96*se2):.3f}-{math.exp(b2+1.96*se2):.3f}",
            f"{p2:.3f}"])
print(f"    OR {math.exp(b2):.3f} "
      f"({math.exp(b2-1.96*se2):.3f}-{math.exp(b2+1.96*se2):.3f}), P = {p2:.3f}")

# -- 3. GEE, exchangeable within species ------------------------------------
dg = df.copy(); dg["grp"] = pd.Categorical(dg.species).codes
m3 = sm.GEE(dg.y, sm.add_constant(dg[["gh"]].astype(float)), groups=dg.grp,
            family=sm.families.Binomial(),
            cov_struct=sm.cov_struct.Exchangeable()).fit()
b3, se3, p3 = m3.params["gh"], m3.bse["gh"], m3.pvalues["gh"]
res.append(["GEE, exchangeable within species", "genome", len(dg),
            f"OR {math.exp(b3):.3f}",
            f"{math.exp(b3-1.96*se3):.3f}-{math.exp(b3+1.96*se3):.3f}",
            f"{p3:.3f}"])
print(f"[3] GEE: OR {math.exp(b3):.3f} "
      f"({math.exp(b3-1.96*se3):.3f}-{math.exp(b3+1.96*se3):.3f}), P = {p3:.3f}")

# -- 4. variance partition (between-species share) --------------------------
grand = df.gh.mean()
ssb = sum(len(v) * (v.gh.mean() - grand) ** 2 for _, v in df.groupby("species"))
sst = ((df.gh - grand) ** 2).sum()
res.append(["Variance partition, GH family count", "species", df.species.nunique(),
            f"{100*ssb/sst:.1f}% between species", "-", "-"])
print(f"[4] between-species variance = {100*ssb/sst:.1f}%")

# -- phylogenetic analyses --------------------------------------------------------
def read_newick(path):
    s = open(path).read().strip().rstrip(";")
    pos = [0]
    def parse():
        node = {"name": "", "len": 0.0, "kids": []}
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                node["kids"].append(parse())
                if s[pos[0]] == ",": pos[0] += 1
                else: break
            pos[0] += 1                      # ')'
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",():": pos[0] += 1
        node["name"] = s[start:pos[0]].strip().strip("'\"")
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1; st = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",()": pos[0] += 1
            node["len"] = float(s[st:pos[0]])
        return node
    return parse()

def tips(n, acc=None):
    acc = [] if acc is None else acc
    if not n["kids"]: acc.append(n)
    else:
        for k in n["kids"]: tips(k, acc)
    return acc

def root_to_tip(n, d=0.0, out=None):
    out = {} if out is None else out
    d += n["len"]
    if not n["kids"]: out[n["name"]] = d
    else:
        for k in n["kids"]: root_to_tip(k, d, out)
    return out

def mrca_depth(n, a, b):
    """Depth to the MRCA of two taxa (length of the shared path)."""
    def path(node, target, d=0.0, acc=None):
        acc = [] if acc is None else acc
        d += node["len"]
        if not node["kids"]:
            return (d, list(acc)) if node["name"] == target else None
        for k in node["kids"]:
            r = path(k, target, d, acc + [node])
            if r: return r
        return None
    pa, pb = path(n, a), path(n, b)
    if not pa or not pb: return 0.0
    shared = 0.0; na, nb = pa[1], pb[1]
    for x, y in zip(na, nb):
        if x is y: shared += x["len"]
        else: break
    return shared

def reroot_outgroup(tree, keyword="Bifidobacterium"):
    """Root at the midpoint of the outgroup (Bifidobacterium) MRCA branch.
    FastTree output is unrooted and the VCV matrix is defined by distance to
    the root, so where the root sits changes the result."""
    def collect(n, path, out):
        if not n["kids"]:
            if keyword.lower() in n["name"].lower(): out.append(list(path) + [n])
        else:
            for k in n["kids"]: collect(k, path + [n], out)
    paths = []; collect(tree, [], paths)
    if len(paths) < 2: return tree
    mrca = tree
    for nodes in zip(*paths):
        if all(x is nodes[0] for x in nodes): mrca = nodes[0]
        else: break
    if mrca is tree: return tree
    # detach the MRCA and place a new root at the midpoint of its branch
    def find_parent(n, target):
        for k in n["kids"]:
            if k is target: return n
            r = find_parent(k, target)
            if r: return r
        return None
    par = find_parent(tree, mrca)
    if par is None: return tree
    half = mrca["len"] / 2.0
    par["kids"] = [k for k in par["kids"] if k is not mrca]
    rest = par if len(par["kids"]) != 1 else par["kids"][0]
    if rest is not par: rest["len"] += par["len"]
    mrca2 = dict(mrca); mrca2["len"] = half
    rest2 = dict(rest); rest2["len"] = rest2["len"] + half
    return {"name": "root", "len": 0.0, "kids": [mrca2, rest2]}


def vcv(tree, names):
    V = np.zeros((len(names), len(names)))
    rt = root_to_tip(tree)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            V[i, j] = rt[a] if i == j else mrca_depth(tree, a, b)
    return V

def pagel_lambda(V, y):
    n = len(y); off = ~np.eye(n, dtype=bool)
    def nll(lam):
        Vl = V.copy(); Vl[off] *= lam
        try: Vi = np.linalg.inv(Vl)
        except np.linalg.LinAlgError: return 1e9
        one = np.ones((n, 1)); yv = y.reshape(-1, 1)
        beta = ((one.T @ Vi @ yv) / (one.T @ Vi @ one)).item()
        r = yv - beta
        s2 = (r.T @ Vi @ r).item() / n
        if s2 <= 0: return 1e9
        sign, ld = np.linalg.slogdet(Vl)
        return 0.5 * (n * math.log(2 * math.pi * s2) + ld + n)
    o = minimize_scalar(nll, bounds=(1e-6, 1.0), method="bounded")
    lam = o.x
    lrt = 2 * (nll(1e-6) - nll(lam))
    p = stats.chi2.sf(max(lrt, 0), 1)
    return lam, max(lrt, 0.0), p

def pgls(V, x, y, lam=1.0):
    n = len(y); off = ~np.eye(n, dtype=bool)
    Vl = V.copy(); Vl[off] *= lam
    Vi = np.linalg.inv(Vl)
    X = np.column_stack([np.ones(n), x])
    XtVi = X.T @ Vi
    beta = np.linalg.solve(XtVi @ X, XtVi @ y)
    r = y - X @ beta
    s2 = (r.T @ Vi @ r).item() / (n - 2)
    cov = s2 * np.linalg.inv(XtVi @ X)
    se = math.sqrt(cov[1, 1])
    t = beta[1] / se
    return beta[1], se, stats.t.sf(abs(t), n - 2) * 2

sp = df.groupby("species").agg(gh_med=("gh", "median"), k=("y", "sum"),
                               n=("y", "count")).reset_index()
sp["prev"] = sp.k / sp.n
sp["asin"] = np.arcsin(np.sqrt(sp.prev))
sp["bin"] = (sp.k > 0).astype(float)

def norm(s): return s.replace(" ", "_").replace("'", "")
for label, path in [("rpoB", "s10_rpoB_species_tree.nwk"),
                    ("16S", os.path.join("twogate", "s16.nwk"))]:
    fp = os.path.join(SRC, path)
    if not os.path.exists(fp):
        print(f"[phylogeny] {label}: tree not found"); continue
    tr = read_newick(fp)
    tn = [t["name"] for t in tips(tr)]
    m = {}
    for s in sp.species:
        for t in tn:
            if norm(s).lower() == norm(t).lower() or norm(s).lower() in norm(t).lower() \
               or norm(t).lower() in norm(s).lower():
                m[s] = t; break
    use = sp[sp.species.isin(m)].copy()
    print(f"\n[phylogeny {label}] tree taxa {len(tn)}, matched {len(use)}/{len(sp)}")
    if len(use) < 5: continue
    names = [m[s] for s in use.species]
    V = vcv(tr, names)
    for trait, col in [("GH family median", "gh_med"),
                       ("tyrDC prevalence (arcsine)", "asin"),
                       ("tyrDC carriage (binary)", "bin")]:
        lam, lrt, p = pagel_lambda(V, use[col].to_numpy(float))
        res.append([f"Pagel lambda, {trait}", "species", len(use),
                    f"lambda {lam:.3f} ({label})", f"LRT {lrt:.2f}", f"{p:.2g}"])
        print(f"   Pagel λ {trait:28s} = {lam:.3f}  LRT {lrt:.2f}  P = {p:.2g}")
    for trait, col in [("prevalence (arcsine)", "asin"), ("carriage (binary)", "bin")]:
        sl, se_, p = pgls(V, use.gh_med.to_numpy(float), use[col].to_numpy(float))
        res.append([f"PGLS, {trait} ~ GH median", "species", len(use),
                    f"slope {sl:+.4f} ({label})",
                    f"{sl-1.96*se_:+.4f}-{sl+1.96*se_:+.4f}", f"{p:.3f}"])
        print(f"   PGLS  {trait:28s} slope {sl:+.4f}  P = {p:.3f}")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["Model or test", "Unit", "n", "Estimate", "95% CI / LRT", "P"])
    w.writerows(res)
print(f"\nsaved {OUT}  ({len(res)}rows)")
