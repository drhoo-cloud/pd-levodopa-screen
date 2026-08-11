#!/usr/bin/env python3
# ============================================================
# 04_validate_controls.py  —  대조군 검증  ★ 여기서 멈출 수 있습니다
# ============================================================
# 실행:  python 04_validate_controls.py
#
# 왜 이 단계가 가장 중요한가
#   전체를 돌려서 "tyrDC 가 거의 없다" 는 결과가 나왔다고 합시다.
#   그것이 정말 없는 것인지, 검색이 작동하지 않은 것인지 구분할 방법이 없습니다.
#   기대 답을 아는 균주 몇 개를 먼저 돌려, 그 답이 나오는지 확인해야 합니다.
#
#   양성 대조가 통과하지 못하면 이 스크립트는 여기서 멈춥니다.
#   참조서열이나 임계값을 고친 뒤 다시 돌리십시오.
# ============================================================

import os
import sys
import csv
import json
import subprocess
import datetime

REF = "refs/gate2_reference.faa"
PROT = "proteins/all_proteins.faa"
OUT = "search/controls.tsv"
LOG = "logs/run_log.txt"

# ------------------------------------------------------------
# 대조군 — accession 은 02 단계 결과에서 확정한 값으로 채우십시오
#   expect: "present"  양성 대조 (온전한 tyrDC 를 가짐)
#           "truncated" 유전자는 있으나 절단 — 판정 세분화 검증용
#           "absent"   음성 대조
# ------------------------------------------------------------
CONTROLS = [
    ("GCF_009734005.1", "Enterococcus faecium SRR24", "present"),
    ("GCF_029023785.1", "Enterococcus faecium DSM 20477", "present"),
    ("GCF_001720945.1", "Enterococcus faecium ISMMS_VRE_1", "present"),
    ("GCF_056485875.1", "Enterococcus faecium s47e1", "present"),
    ("GCF_056485845.1", "Enterococcus faecium s47v1", "present"),
    ("GCF_009697285.1", "Enterococcus faecium VRE", "present"),
    ("GCF_019977575.1", "Enterococcus faecium AA622", "present"),
    ("GCF_900639535.1", "Enterococcus faecium -", "present"),
    ("GCF_056485905.1", "Enterococcus faecium p463s", "present"),
    ("GCF_002007625.1", "Enterococcus faecium 2014-VREF-41", "present"),
    ("GCF_056485925.1", "Enterococcus faecium p344", "present"),
    ("GCF_900639715.1", "Enterococcus faecium -", "present"),
    ("GCF_029024925.1", "Enterococcus faecalis DSM 20478", "present"),
    ("GCF_022869705.1", "Enterococcus faecalis PartL-Efaecalis-RM8376", "present"),
    ("GCF_021610105.1", "Enterococcus faecalis UK045", "present"),
    ("GCF_050485625.1", "Enterococcus faecalis C2198", "present"),
    ("GCF_018986755.2", "Enterococcus faecalis 1207/14", "present"),
    ("GCF_050485045.1", "Enterococcus faecalis C8816", "present"),
    ("GCF_050485635.1", "Enterococcus faecalis C2196", "present"),
    ("GCF_050485055.1", "Enterococcus faecalis C8811", "present"),
    ("GCF_006494835.1", "Enterococcus faecalis VE14089", "present"),
    ("GCF_028131645.1", "Enterococcus faecalis SVR2330", "present"),
    ("GCF_006494855.1", "Enterococcus faecalis VE18379", "present"),
    ("GCF_055814325.1", "Enterococcus faecalis 1604D004", "present"),
    ("GCF_000007785.1", "Enterococcus faecalis V583", "present"),
    ("GCF_000203855.3", "Lactiplantibacillus plantarum WCFS1", "absent"),
    ("GCF_001005805.1", "Lactiplantibacillus plantarum PS128", "absent"),
]# 임계값 — 사전 확정. 바꾸면 로그에 남기고 이유를 적으십시오.
ID_PRESENT, COV_PRESENT = 60.0, 80.0
ID_AMBIG, COV_AMBIG = 40.0, 50.0
MARGIN = 10.0          # 최상위 계열과 차상위 계열의 점수 차이 최소값


def run_diamond():
    os.makedirs("search", exist_ok=True)
    if not os.path.exists("search/ref.dmnd"):
        print("  DIAMOND 데이터베이스 만드는 중")
        subprocess.run(["diamond", "makedb", "--in", REF, "-d", "search/ref"],
                       check=True, capture_output=True)
    print("  검색 중 (대조군)")
    subprocess.run([
        "diamond", "blastp", "-q", PROT, "-d", "search/ref", "-o", OUT,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "qlen", "slen", "evalue", "bitscore",
        "--evalue", "1e-5", "--max-target-seqs", "25", "--quiet",
    ], check=True)


CALL_FAMILIES = ('TARGET_tyrDC', 'TARGET_aadc')


def call_family(hits):
    """한 단백질의 히트들을 보고 어느 계열인지, 어떤 판정인지 정합니다"""
    best = {}
    for h in hits:
        fam = h["sseqid"].split("|")[0]          # TARGET_tyrDC / DECOY_gadB ...
        if fam not in best or h["bitscore"] > best[fam]["bitscore"]:
            best[fam] = h
    ranked = sorted(best.items(), key=lambda kv: -kv[1]["bitscore"])
    if not ranked:
        return "absent", None, 0.0
    top_fam, top = ranked[0]
    second = ranked[1][1]["bitscore"] if len(ranked) > 1 else 0.0
    margin = top["bitscore"] - second

    if top_fam not in CALL_FAMILIES:
        return "absent", top_fam, margin        # 최상위가 decoy → 표적 아님
    cov = 100.0 * top["length"] / max(top["slen"], 1)
    if top["pident"] >= ID_PRESENT and cov >= COV_PRESENT and margin >= MARGIN:
        return "present", top_fam, margin
    if top["pident"] >= ID_AMBIG and cov >= COV_AMBIG:
        return "ambiguous", top_fam, margin
    return "absent", top_fam, margin


def main():
    for f in (REF, PROT):
        if not os.path.exists(f):
            sys.exit(f"파일이 없습니다: {f}  (앞 단계를 먼저 실행하십시오)")

    todo = [c for c in CONTROLS if c[0]]
    if len(todo) < 3:
        print("★ 대조군 accession 이 3건 미만입니다.")
        print("  04_validate_controls.py 의 CONTROLS 목록을 채운 뒤 다시 실행하십시오.")
        print("  panel/strain_resolution.tsv 에서 확정된 accession 을 찾을 수 있습니다.")
        sys.exit(1)

    print("=== 4단계: 대조군 검증 ===")
    run_diamond()

    # accession 별로 히트 모으기
    by_acc = {}
    with open(OUT) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            acc = p[0].split("|")[0]
            by_acc.setdefault(acc, []).append({
                "qseqid": p[0], "sseqid": p[1], "pident": float(p[2]),
                "length": int(p[3]), "qlen": int(p[4]), "slen": int(p[5]),
                "evalue": float(p[6]), "bitscore": float(p[7])})

    print("\n  결과")
    print("  " + "-" * 66)
    ok, fail = 0, []
    for acc, name, expect in todo:
        hits = by_acc.get(acc, [])
        # 이 유전체의 단백질을 하나씩 판정해 가장 강한 결론을 취합니다
        calls = []
        by_prot = {}
        for h in hits:
            by_prot.setdefault(h["qseqid"], []).append(h)
        for _, hs in by_prot.items():
            calls.append(call_family(hs)[0])
        if "present" in calls:
            got = "present"
        elif "ambiguous" in calls:
            got = "ambiguous"
        else:
            got = "absent"

        # truncated 기대는 present 또는 ambiguous 로 잡히면 통과로 봅니다
        passed = (got == expect) or (expect == "truncated" and got in ("present", "ambiguous"))
        mark = "OK  " if passed else "실패"
        print(f"  {mark} {name:36s} 기대={expect:9s} 실제={got}")
        if passed:
            ok += 1
        else:
            fail.append((name, expect, got))

    print("  " + "-" * 66)
    with open(LOG, "a") as lg:
        lg.write(f"\n[04_validate_controls] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  대조군 {len(todo)} · 통과 {ok} · 실패 {len(fail)}\n")
        lg.write(f"  임계값 present id>={ID_PRESENT} cov>={COV_PRESENT} margin>={MARGIN}\n")
        for n, e, g in fail:
            lg.write(f"  실패: {n} 기대={e} 실제={g}\n")

    if fail:
        print("\n" + "=" * 66)
        print("★ 여기서 멈춥니다. 전체 분석으로 넘어가지 마십시오.")
        print("=" * 66)
        print("  확인할 것:")
        print("   1) refs/gate2_reference.faa 에 해당 계열 서열이 실제로 들어 있는가")
        print("   2) 임계값이 너무 엄격한가 (ID_PRESENT, COV_PRESENT, MARGIN)")
        print("   3) 그 유전체의 protein.faa 가 제대로 받아졌는가")
        print("  고친 뒤 이 스크립트를 다시 실행하십시오.")
        sys.exit(2)

    print("\n=== 4단계 통과 ===")
    print("  대조군이 모두 기대한 판정을 반환했습니다.")
    print("다음:  bash 05_search_all.sh")


if __name__ == "__main__":
    main()
