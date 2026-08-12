#!/usr/bin/env python3
# ============================================================
# 06_call_gates.py  —  유전체별 판정
# ============================================================
# 실행:  python 06_call_gates.py
#
# 무엇을 하나
#   검색 결과를 읽어 유전체마다 다음을 정합니다.
#     Gate 2 : present / ambiguous / absent
#              present 는 다시 intact / truncated / point-variant 로 나눕니다
#     Gate 1 : GH family 개수와 pass/fail
#
# ★ 반드시 지키는 규칙
#   ambiguous 를 absent 에 합치지 않습니다.
#   합치는 순간 "작동하는 스크린" 이 "안심시키는 스크린" 으로 바뀝니다.
# ============================================================

import os
import csv
import json
import datetime
from collections import defaultdict

HITS = "search/gate2_hits.tsv"
DBCAN = "search/dbcan/overview.txt"
OUT = "results/per_genome_calls.tsv"
LOG = "logs/run_log.txt"

# 사전 확정 임계값 — 04 단계와 같은 값이어야 합니다
ID_PRESENT, COV_PRESENT = 60.0, 80.0
ID_AMBIG, COV_AMBIG = 40.0, 50.0
MARGIN = 10.0
GH_MIN = 15                     # Gate 1 통과 최소 GH family 수
TRUNC_RATIO = 0.80              # 참조 길이의 80% 미만이면 절단으로 봅니다


def load_hits(path):
    by_prot = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            by_prot[p[0]].append({
                "sseqid": p[1], "pident": float(p[2]), "length": int(p[3]),
                "qlen": int(p[4]), "slen": int(p[5]),
                "evalue": float(p[6]), "bitscore": float(p[7])})
    return by_prot


CALL_FAMILIES = ('TARGET_tyrDC', 'TARGET_aadc')


def classify(hits):
    """한 단백질에 대해 (판정, 계열, 변이유형, 마진) 반환"""
    best = {}
    for h in hits:
        fam = h["sseqid"].split("|")[0]
        if fam not in best or h["bitscore"] > best[fam]["bitscore"]:
            best[fam] = h
    ranked = sorted(best.items(), key=lambda kv: -kv[1]["bitscore"])
    if not ranked:
        return "absent", "", "", 0.0
    fam, top = ranked[0]
    second = ranked[1][1]["bitscore"] if len(ranked) > 1 else 0.0
    margin = round(top["bitscore"] - second, 1)

    # 최상위가 decoy 면 표적이 아닙니다 — 이것이 오배정을 막는 지점입니다
    if fam not in CALL_FAMILIES:
        return "absent", fam, "decoy_top", margin

    cov = 100.0 * top["length"] / max(top["slen"], 1)
    len_ratio = top["qlen"] / max(top["slen"], 1)

    if top["pident"] >= ID_PRESENT and cov >= COV_PRESENT and margin >= MARGIN:
        if len_ratio < TRUNC_RATIO:
            variant = "truncated"          # 질의 단백질이 참조보다 뚜렷하게 짧음
        elif top["pident"] < 95.0:
            variant = "point-variant"      # 온전하지만 치환이 많음 → 활성 확인 필요
        else:
            variant = "intact"
        return "present", fam, variant, margin

    if top["pident"] >= ID_AMBIG and cov >= COV_AMBIG:
        return "ambiguous", fam, "", margin

    return "absent", fam, "", margin


def load_gh_counts(path):
    """dbCAN overview 에서 유전체별 GH family 종류 수를 셉니다"""
    counts = defaultdict(set)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            acc = p[0].split("|")[0]
            for field in p[1:]:
                for tok in str(field).replace("+", " ").replace(";", " ").split():
                    if tok.startswith("GH"):
                        counts[acc].add(tok.split("_")[0].split("(")[0])
    return {k: len(v) for k, v in counts.items()}


def main():
    os.makedirs("results", exist_ok=True)
    if not os.path.exists(HITS):
        raise SystemExit(f"파일이 없습니다: {HITS}  (05 단계를 먼저 실행하십시오)")

    by_prot = load_hits(HITS)
    # 히트가 없는 유전체도 분모에 포함한다 (명백한 음성)
    ALL_ACC = set()
    import io
    with open("proteins/all_proteins.faa") as _f:
        for _l in _f:
            if _l.startswith(">"):
                ALL_ACC.add(_l[1:].split("|")[0])
    print(f"  전체 유전체 {len(ALL_ACC)}건 (히트 유무 무관)")
    gh = load_gh_counts(DBCAN)

    # 유전체 단위로 모읍니다
    genome = defaultdict(lambda: {"present": [], "ambiguous": [], "absent": 0})
    for qseqid, hits in by_prot.items():
        acc = qseqid.split("|")[0]
        call, fam, variant, margin = classify(hits)
        if call == "present":
            genome[acc]["present"].append((qseqid, fam, variant, margin))
        elif call == "ambiguous":
            genome[acc]["ambiguous"].append((qseqid, fam, margin))
        else:
            genome[acc]["absent"] += 1

    rows, tally = [], defaultdict(int)
    for acc in sorted(ALL_ACC):
        _ = genome[acc]
        g = genome[acc]
        if g["present"]:
            g2 = "present"
            variants = ";".join(sorted({v for _, _, v, _ in g["present"]}))
            loci = ";".join(q for q, _, _, _ in g["present"])
        elif g["ambiguous"]:
            g2, variants = "ambiguous", ""
            loci = ";".join(q for q, _, _ in g["ambiguous"])
        else:
            g2, variants, loci = "absent", "", ""

        n_gh = gh.get(acc, "")
        g1 = "" if n_gh == "" else ("pass" if n_gh >= GH_MIN else "fail")

        tally[g2] += 1
        rows.append({
            "assembly_accession": acc,
            "gate2_call": g2,
            "gate2_variant": variants,
            "gate2_loci": loci,
            "gate1_gh_families": n_gh,
            "gate1_call": g1,
        })

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print("=== 6단계: 판정 결과 ===")
    print(f"  유전체 {len(rows)} 건")
    for k in ("present", "ambiguous", "absent"):
        print(f"    Gate 2 {k:10s} {tally[k]:5d}")
    if gh:
        n_pass = sum(1 for r in rows if r["gate1_call"] == "pass")
        print(f"    Gate 1 pass       {n_pass:5d}  (GH >= {GH_MIN})")
    else:
        print("    Gate 1 미실행 — dbCAN 을 설치하고 05 단계를 다시 돌리십시오")

    with open(LOG, "a") as lg:
        lg.write(f"\n[06_call_gates] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  유전체 {len(rows)} · present {tally['present']} "
                 f"· ambiguous {tally['ambiguous']} · absent {tally['absent']}\n")
        lg.write(f"  임계값 id>={ID_PRESENT}/{ID_AMBIG} cov>={COV_PRESENT}/{COV_AMBIG} "
                 f"margin>={MARGIN} GH_MIN={GH_MIN}\n")

    print(f"\n  → {OUT}")
    print("다음:  python 07_summarize.py")


if __name__ == "__main__":
    main()
