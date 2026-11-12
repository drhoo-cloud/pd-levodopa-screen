#!/usr/bin/env python3
"""Recompute every headline number of the paper from the deposited files.

Usage:   python3 verify_all.py          (run from the repository root)

No third-party packages are required. Each line prints the value stated in
the paper, the value recomputed here, and whether they agree.
"""
import csv, math, os, sys, collections, statistics as st, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
def rd(p):
    with open(os.path.join(HERE, p), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

OK = FAIL = 0
def chk(label, stated, got, ok=None):
    global OK, FAIL
    good = (str(stated) == str(got)) if ok is None else ok
    OK, FAIL = OK + good, FAIL + (not good)
    print(f"  [{'OK ' if good else 'DIFF'}] {label:<52} paper {str(stated):>14}   file {str(got):>14}")

def wilson_upper(n, z=1.959963985):
    d = 1 + z * z / n
    return ((z * z / (2 * n)) / d + (z / d) * math.sqrt(z * z / (4 * n * n))) * 100

print("=" * 104)
print("Verification of the main-text numbers against the deposited files")
print("=" * 104)

# ---- Gate 2 census -------------------------------------------------------
print("\nSection 3  Gate 2 census")
g2 = [r for r in rd("data/gate2_calls_per_genome.tsv") if r["in_census"] == "yes"]
chk("genomes screened at Gate 2", "19771", f"{len(g2)}")
pos = [r for r in g2 if r["gate2_call"].lower() == "present"]
chk("genomes positive", "11365", f"{len(pos)}")
chk("pooled prevalence, %", "57.5", f"{100*len(pos)/len(g2):.1f}")
sp = collections.defaultdict(lambda: [0, 0])
for r in g2:
    sp[r["species"]][1] += 1
    if r["gate2_call"].lower() == "present": sp[r["species"]][0] += 1
chk("species with the locus", "6", sum(1 for k, n in sp.values() if k))
chk("species with no detection", "16", sum(1 for k, n in sp.values() if not k))
ctx = rd("data/locus_context_by_species.tsv")
chk("species with the canonical operon", "5", sum(1 for r in ctx if int(r["canonical"]) > 0))
lb = sp["Levilactobacillus brevis"]
chk("L. brevis prevalence, %", "80.1", f"{100*lb[0]/lb[1]:.1f}")
lp = sp["Lactiplantibacillus plantarum"]
chk("L. plantarum carriers / genomes", "0 / 1665", f"{lp[0]} / {lp[1]}")
zb = {s: wilson_upper(n) for s, (k, n) in sp.items() if k == 0}
chk("lowest Wilson upper bound, %", "0.23", f"{min(zb.values()):.2f}")
chk("highest Wilson upper bound, %", "9.89", f"{max(zb.values()):.2f}")
pp = [r for r in ctx if "pentosaceus" in r["species"]][0]
chk("P. pentosaceus calls without the operon", "5", pp["no_operon"])

# ---- Gate 1 --------------------------------------------------------------
print("\nSection 3  Gate 1 census")
g1 = rd("data/gate1_gh_families_per_genome.tsv")
chk("genomes screened at Gate 1", "19765", f"{len(g1)}")
fam = collections.defaultdict(list)
for r in g1: fam[r["species"]].append(int(r["gh_families"]))
med = {s: st.median(v) for s, v in fam.items()}
chk("lowest species median, GH families", "7", f"{min(med.values()):.0f}")
chk("highest species median, GH families", "28", f"{max(med.values()):.0f}")
def q(v, p):
    v = sorted(v); k = (len(v) - 1) * p; f = int(k)
    return v[f] + (v[min(f + 1, len(v) - 1)] - v[f]) * (k - f)
iqr = {s: (q(v, .25), q(v, .75)) for s, v in fam.items()}
pr = list(itertools.combinations(sorted(fam), 2))
nov = sum(1 for a, b in pr if iqr[a][1] < iqr[b][0] or iqr[b][1] < iqr[a][0])
chk("species pairs with non-overlapping IQR, %", "81", f"{100*nov/len(pr):.0f}")
caz = rd("data/cazyme_all_categories.tsv")
chk("species lacking polysaccharide lyases", "20",
    sum(1 for r in caz if int(float(r["PL_families"])) == 0))

# ---- L. plantarum -------------------------------------------------------
print("\nSection 3  Absence in L. plantarum")
ph = rd("data/lp_prophage_census.tsv")
chk("genomes with a prophage, %", "98.3",
    f"{100*sum(1 for r in ph if int(r['prophage_regions'] or 0) > 0)/len(ph):.1f}")
chk("prophage regions, total", "7305",
    sum(int(r["prophage_regions"] or 0) for r in ph))
pf = rd("data/lp_pf00282_in_prophage.tsv")
chk("genomes carrying the PF00282 family, %", "98.62", f"{100*len(pf)/len(ph):.2f}")
chk("PF00282 proteins, total", "1654",
    sum(int(r["pf00282_proteins"] or 0) for r in pf))
chk("PF00282 proteins inside a prophage", "27",
    sum(int(r["pf00282_in_prophage"] or 0) for r in pf))

# ---- dereplication ------------------------------------------------------
print("\nSection 4  Dereplication")
dr = rd("data/lp_dereplication_census.tsv")[0]
chk("genomes dereplicated", "1665", dr["n_genomes"])
chk("independent lineages", "161", dr["n_clusters"])
chk("Wilson upper bound after dereplication, %", "2.33",
    f"{wilson_upper(int(dr['n_clusters'])):.2f}")
mem = rd("data/lp_cluster_membership.tsv")
chk("cluster rows match the lineage count", dr["n_clusters"], f"{len(mem)}")
chk("members sum to the genome count", dr["n_genomes"],
    f"{sum(int(r['n_members']) for r in mem)}")

# ---- isolation source ---------------------------------------------------
iso = rd("data/isolation_source.tsv")
grp = collections.defaultdict(lambda: [0, 0])
for r in iso:
    grp[r["src_cat"]][1] += 1
    if r["tyrdc"] not in ("0", ""): grp[r["src_cat"]][0] += 1
a, n = grp["probiotic / product"]
chk("prevalence among probiotic isolates, %", "1.28", f"{100*a/n:.2f}")
a, n = grp["environment / industrial"]
chk("prevalence among environmental isolates, %", "34.53", f"{100*a/n:.2f}")

# ---- profile validation -------------------------------------------------
print("\nSection 4  What the profile does not see")
sv = rd("data/profile_validation_scores.tsv")
bac = [r for r in sv if "Bacill" in r["phylum"]]
out = [r for r in sv if "Bacill" not in r["phylum"]]
chk("reviewed sequences within Bacillota", "147", f"{len(bac)}")
chk("of which detected", "145", sum(1 for r in bac if r["call"] == "present"))
chk("reviewed sequences outside Bacillota", "211", f"{len(out)}")
chk("median identity outside Bacillota, %", "32.7",
    f"{st.median(float(r['max_id_to_ref']) for r in out):.1f}")
det = [r for r in out if r["call"] == "present"]
chk("of which clear the threshold", "4", f"{len(det)}")
chk("lowest score among those, %", "59.5",
    f"{min(float(r['score_pct']) for r in det):.1f}")
chk("highest score among those, %", "74.3",
    f"{max(float(r['score_pct']) for r in det):.1f}")
nd = max(out, key=lambda r: float(r["score_pct"]) if r["call"] != "present" else -1)
chk("highest excluded score, %", "50.00", f"{float(nd['score_pct']):.2f}")
bm = rd("data/profile_benchmark_scores.tsv")
neg = [r for r in bm if r["truth"] == "neg"]
chk("related decarboxylases tested", "166", f"{len(neg)}")
chk("of which called positive", "0", sum(1 for r in neg if float(r["score_pct"]) > 50))
chk("highest negative score, %", "4.72",
    f"{max(float(r['score_pct']) for r in neg):.2f}")
lowdet = min(float(r["score_pct"]) for r in bac if r["call"] == "present")
chk("threshold-invariant interval, points", "53.7",
    f"{lowdet - max(float(r['score_pct']) for r in neg):.1f}")

# ---- pre-gate -----------------------------------------------------------
print("\nSection 2  Pre-gate")
pg = rd("data/pregate_per_genome.tsv")
chk("control genomes", "25", f"{len(pg)}")
chk("failing the pre-gate", "25", sum(1 for r in pg if r["pre_gate_call"].startswith("FAIL")))
pc = rd("data/pregate_context.tsv")
chk("determinant hits in total", "822", f"{len(pc)}")
van = [x for x in pc if "van" in x["gene_product"].lower()]
oth = [x for x in pc if "van" not in x["gene_product"].lower()]
chk("van operon gene hits", "123", f"{len(van)}")
for w, lab in (("10kb", "10 kb"), ("5kb", "5 kb")):
    c = f"mobile_adjacent_{w}"
    a = sum(1 for x in van if x[c] == "1")
    o = sum(1 for x in oth if x[c] == "1")
    stated = {"10kb": ("100.0", "45.8"), "5kb": ("90.2", "35.1")}[w]
    chk(f"van mobile-adjacent within {lab}, %", stated[0], f"{100*a/len(van):.1f}")
    chk(f"others mobile-adjacent within {lab}, %", stated[1], f"{100*o/len(oth):.1f}")

# ---- two gates ----------------------------------------------------------
print("\nTable S4b  The two gates are not interchangeable")
tg = {r["Model or test"]: r for r in rd("data/twogate_census_stats.tsv")}
chk("logistic without a species term, n", "19765", tg["Logistic, no species term"]["n"])
chk("odds ratio per GH family", "OR 1.297 per GH family",
    tg["Logistic, no species term"]["Estimate"])
chk("species fixed effects, n", "12273", tg["Logistic, species fixed effects"]["n"])
chk("species fixed effects, P", "0.107", tg["Logistic, species fixed effects"]["P"])
chk("variance between species", "88.6% between species",
    tg["Variance partition, GH family count"]["Estimate"])

print("\n" + "=" * 104)
print(f"  {OK} reproduced, {FAIL} differing")
print("""
One number is not recomputed here. The eleven-trial disclosure audit of
Section 1 and Table S1 is a manual reading of the publications, tabulated in
the Supplementary Information; data/disclosure_audit_named_strains.tsv holds
the 22 components that are resolved to strain.""")
sys.exit(1 if FAIL else 0)
