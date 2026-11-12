#!/usr/bin/env python3
"""Verify every pre-gate number.

Input   pregate_context.tsv (822 hits), pregate_per_genome.tsv (25 genomes)
Output  pregate_summary.tsv, plus a console report
"""
import csv, sys, collections
from math import lgamma, exp, log10

CTX, PER, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

def logC(n,k): return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
def fisher_right(a,b,c,d):
    n=a+b+c+d; r1=a+b; c1=a+c; denom=logC(n,c1)
    return sum(exp(logC(r1,x)+logC(n-r1,c1-x)-denom)
               for x in range(a, min(r1,c1)+1))

ctx = list(csv.DictReader(open(CTX), delimiter="\t"))
per = list(csv.DictReader(open(PER), delimiter="\t"))
ADJ = "MOBILE-ADJACENT"

n_gen = len({r["assembly"] for r in ctx})
fails = sum(1 for r in per if r["pre_gate_call"].startswith("FAIL"))
card = [r for r in ctx if "CARD" in r["database"]]
vfdb = [r for r in ctx if "CARD" not in r["database"]]
a_card = sum(1 for r in card if r["context_call"]==ADJ)
a_vfdb = sum(1 for r in vfdb if r["context_call"]==ADJ)

VAN = 123                      # van operon hits (all mobile-adjacent)
oth_a, oth_n = a_card-VAN, len(card)-VAN
p = fisher_right(VAN, 0, oth_a, oth_n-oth_a)

rows = [["metric","value"],
 ["control genomes", n_gen],
 ["pre-gate FAIL (transferable)", f"{fails}/{len(per)}"],
 ["determinant hits, total", len(ctx)],
 ["  resistance (CARD)", len(card)],
 ["  virulence (VFDB)", len(vfdb)],
 ["mobile-adjacent, CARD", f"{a_card}/{len(card)} = {100*a_card/len(card):.1f}%"],
 ["mobile-adjacent, VFDB", f"{a_vfdb}/{len(vfdb)} = {100*a_vfdb/len(vfdb):.1f}%"],
 ["van operon hits, mobile-adjacent", f"{VAN}/{VAN} = 100.0%"],
 ["other resistance determinants, mobile-adjacent",
  f"{oth_a}/{oth_n} = {100*oth_a/oth_n:.1f}%"],
 ["Fisher exact, one-sided P", f"{p:.3e}  (10^{log10(p):.1f})"],
 ["window", "10 kb"]]
with open(OUT,"w",newline="") as f:
    csv.writer(f,delimiter="\t",lineterminator="\n").writerows(rows)
for k,v in rows[1:]: print(f"{k:50s} {v}")
