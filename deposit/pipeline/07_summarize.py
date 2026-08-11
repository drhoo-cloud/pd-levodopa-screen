#!/usr/bin/env python3
# ============================================================
# 07_summarize.py  —  집계와 원고용 표 만들기
# ============================================================
# 실행:  python 07_summarize.py
#
# 무엇을 하나
#   유전체별 판정을 층(niche, 종, 패널 그룹)별로 집계하고
#   Wilson 95% 신뢰구간과 보수적 상한을 함께 계산해 표로 만듭니다.
#
# ★ 왜 Wilson 인가
#   0/450 을 "0%" 로 쓰면 안 됩니다. 표본이 450건일 뿐 없다는 증거가 아닙니다.
#   Wilson 상한을 함께 써야 "0.82% 이하" 라는 정확한 진술이 됩니다.
#
# ★ 보수적 상한
#   ambiguous 를 전부 양성으로 세었을 때의 값도 같이 냅니다.
#   심사자가 반드시 묻는 질문이고, 미리 답해 두면 논거가 강해집니다.
# ============================================================

import os
import csv
import math
import datetime
from collections import defaultdict

CALLS = "results/per_genome_calls.tsv"
META = "panel/assemblies_all.tsv"
STRATA = "panel/strata.tsv"          # 선택: assembly_accession \t stratum
LOG = "logs/run_log.txt"


def wilson(k, n, z=1.96):
    """Wilson 점수 신뢰구간 — k=0 이어도 상한이 나옵니다"""
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * p, 100 * max(0.0, c - h), 100 * min(1.0, c + h))


def fmt(k, n):
    p, lo, hi = wilson(k, n)
    if k == 0:
        return f"0/{n} — 0%, 95% 상한 {hi:.2f}%"
    return f"{k}/{n} — {p:.3f}% [{lo:.3f}–{hi:.3f}%]"


def main():
    if not os.path.exists(CALLS):
        raise SystemExit(f"파일이 없습니다: {CALLS}  (06 단계를 먼저 실행하십시오)")

    calls = {}
    with open(CALLS) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            calls[row["assembly_accession"]] = row

    # 종 정보 붙이기
    species = {}
    if os.path.exists(META):
        with open(META) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = row.get("organism", "")
                species[row["accession"]] = " ".join(name.split()[:2])

    # 층 정보 (선택)
    strata = {}
    if os.path.exists(STRATA):
        with open(STRATA) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    strata[p[0]] = p[1]

    os.makedirs("results", exist_ok=True)

    # ------------------------------------------------------------
    # 표 A — 종별 집계
    # ------------------------------------------------------------
    by_sp = defaultdict(lambda: {"n": 0, "present": 0, "ambiguous": 0, "gh": []})
    for acc, r in calls.items():
        sp = species.get(acc, "(unknown)")
        d = by_sp[sp]
        d["n"] += 1
        if r["gate2_call"] == "present":
            d["present"] += 1
        elif r["gate2_call"] == "ambiguous":
            d["ambiguous"] += 1
        if r["gate1_gh_families"]:
            d["gh"].append(int(r["gate1_gh_families"]))

    rows_a = []
    for sp in sorted(by_sp, key=lambda s: -by_sp[s]["n"]):
        d = by_sp[sp]
        k, a, n = d["present"], d["ambiguous"], d["n"]
        gh = sorted(d["gh"])
        med = gh[len(gh)//2] if gh else ""
        rows_a.append({
            "species": sp, "n": n,
            "tyrDC_present": k,
            "prevalence_wilson": fmt(k, n),
            "ambiguous": a,
            "conservative_upper": fmt(k + a, n),      # ambiguous 를 양성으로 셈
            "GH_median": med,
        })

    out_a = "results/table_by_species.tsv"
    with open(out_a, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_a[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows_a)

    # ------------------------------------------------------------
    # 표 B — 층별 집계 (분리원 등). strata.tsv 가 있을 때만
    # ------------------------------------------------------------
    rows_b = []
    if strata:
        by_st = defaultdict(lambda: {"n": 0, "present": 0, "ambiguous": 0})
        for acc, r in calls.items():
            st = strata.get(acc, "UNKNOWN")
            d = by_st[st]
            d["n"] += 1
            if r["gate2_call"] == "present":
                d["present"] += 1
            elif r["gate2_call"] == "ambiguous":
                d["ambiguous"] += 1
        for st in sorted(by_st, key=lambda s: -by_st[s]["n"]):
            d = by_st[st]
            rows_b.append({
                "stratum": st, "n": d["n"],
                "tyrDC_present": d["present"],
                "prevalence_wilson": fmt(d["present"], d["n"]),
                "ambiguous": d["ambiguous"],
                "conservative_upper": fmt(d["present"] + d["ambiguous"], d["n"]),
            })
        out_b = "results/table_by_stratum.tsv"
        with open(out_b, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_b[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(rows_b)

    # ------------------------------------------------------------
    # 화면 출력
    # ------------------------------------------------------------
    tot = len(calls)
    tp = sum(1 for r in calls.values() if r["gate2_call"] == "present")
    ta = sum(1 for r in calls.values() if r["gate2_call"] == "ambiguous")
    print("=== 7단계: 집계 ===")
    print(f"  전체 {tot} 유전체")
    print(f"  Gate 2 present    {fmt(tp, tot)}")
    print(f"  Gate 2 ambiguous  {ta}")
    print(f"  보수적 상한       {fmt(tp + ta, tot)}")
    print()
    print("  종별 (상위 8)")
    for r in rows_a[:8]:
        print(f"    {r['species'][:34]:36s} {r['prevalence_wilson']}")
    if rows_b:
        print()
        print("  층별")
        for r in rows_b:
            print(f"    {r['stratum'][:34]:36s} {r['prevalence_wilson']}")

    with open(LOG, "a") as lg:
        lg.write(f"\n[07_summarize] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  전체 {tot} · present {tp} · ambiguous {ta}\n")

    print(f"\n  → {out_a}")
    if rows_b:
        print(f"  → results/table_by_stratum.tsv")
    print()
    print("=== 완료 ===")
    print("원고에 넣을 때 지킬 것")
    print("  · 0 건인 층은 '0%' 가 아니라 '0%, 95% 상한 X%' 로 씁니다")
    print("  · ambiguous 는 별도 열로 보고하고 absent 에 합치지 않습니다")
    print("  · logs/run_log.txt 를 Supplementary 에 그대로 첨부합니다")


if __name__ == "__main__":
    main()
