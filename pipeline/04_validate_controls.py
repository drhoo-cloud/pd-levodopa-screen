#!/usr/bin/env python3
# ============================================================
# 04_validate_controls.py  —  대조군 검증  ★ 여기서 멈출 수 있습니다
# ============================================================
# 실행:  python 04_validate_controls.py
#
# 왜 이 단계가 가장 중요한가
#   전체를 돌려서 "tyrDC 가 거의 없다" 는 결과가 나왔다고 합시다.
#   그것이 정말 없는 것인지, 검색이 작동하지 않은 것인지 구분할 방법이 없습니다.
#   기대 답을 아는 균주를 먼저 돌려, 그 답이 나오는지 확인해야 합니다.
#
#   양성 대조가 통과하지 못하면 이 스크립트는 여기서 멈춥니다.
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점
# ------------------------------------------------------------
#   DIAMOND + identity/coverage/margin → hmmsearch + 정규화 점수 하나
#   ambiguous 범주 삭제
#
#   양성 대조 두 건은 참조 세트 안에 들어 있는 균주입니다.
#     Enterococcus faecalis V583   (EF_0634, UniProt Q838D6)
#     Levilactobacillus brevis     (UniProt J7GQ11)
#   둘 다 TyrDC 활성이 실험으로 확인된 것이고, 구조도 풀려 있습니다.
#   이 둘이 100% 근처로 나오지 않으면 프로파일이나 self-score 가 잘못된 것입니다.
# ============================================================

import os
import sys
import csv
import subprocess
import datetime

HMM = "refs/tyrdc.hmm"
SELF_F = "refs/self_score.txt"
PROT = "proteins/all_proteins.faa"
OUT = "search/controls.tsv"
DOMTBL = "search/controls.domtbl"
LOG = "logs/run_log.txt"

# 06 단계와 반드시 같은 값이어야 합니다
CUTOFF = 50.0
EVALUE = "1e-5"

# ------------------------------------------------------------
# 대조군 — accession 은 02 단계 결과에서 확정한 값입니다
#   expect: "present"  양성 대조
#           "absent"   음성 대조
# ------------------------------------------------------------
CONTROLS = [
    ("GCF_000007785.1", "Enterococcus faecalis V583  ★참조세트", "present"),
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
    ("GCF_022869705.1", "Enterococcus faecalis RM8376", "present"),
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
    ("GCF_000203855.3", "Lactiplantibacillus plantarum WCFS1", "absent"),
    ("GCF_001005805.1", "Lactiplantibacillus plantarum PS128", "absent"),
]


def run_hmmsearch():
    os.makedirs("search", exist_ok=True)
    print("  hmmsearch 실행 (대조군)")
    subprocess.run([
        "hmmsearch", "--domtblout", DOMTBL, "-o", "/dev/null",
        "-E", EVALUE, HMM, PROT,
    ], check=True)


def best_by_genome(path, self_score):
    """유전체마다 가장 높은 정규화 점수를 돌려줍니다."""
    best = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 8:
            continue
        acc = p[0].split("|")[0]
        pct = 100.0 * float(p[7]) / self_score
        if pct > best.get(acc, 0.0):
            best[acc] = pct
    return best


def main():
    for f in (HMM, SELF_F, PROT):
        if not os.path.exists(f):
            sys.exit(f"파일이 없습니다: {f}  (앞 단계를 먼저 실행하십시오)")

    self_score = float(open(SELF_F).read().strip())

    todo = [c for c in CONTROLS if c[0]]
    if len(todo) < 3:
        sys.exit("★ 대조군이 3건 미만입니다. CONTROLS 목록을 채우십시오.")

    print("=== 4단계: 대조군 검증 ===")
    print(f"  프로파일 {HMM}")
    print(f"  self-score {self_score} · 임계 {CUTOFF}%\n")
    run_hmmsearch()

    best = best_by_genome(DOMTBL, self_score)

    print("\n  결과")
    print("  " + "-" * 72)
    ok, fail = 0, []
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["assembly_accession", "name", "expected", "score_pct", "observed"])
        for acc, name, expect in todo:
            pct = best.get(acc, 0.0)
            got = "present" if pct > CUTOFF else "absent"
            passed = (got == expect)
            mark = "OK  " if passed else "실패"
            print(f"  {mark} {name:38s} 기대={expect:8s} 점수={pct:6.2f}%  실제={got}")
            w.writerow([acc, name, expect, f"{pct:.2f}", got])
            if passed:
                ok += 1
            else:
                fail.append((name, expect, got, pct))
    print("  " + "-" * 72)

    # 참조 세트에 들어 있는 균주는 100% 근처여야 합니다
    v583 = best.get("GCF_000007785.1", 0.0)
    if v583 and v583 < 95.0:
        print(f"\n  ★ 주의: E. faecalis V583 가 {v583:.2f}% 입니다.")
        print("     이 균주의 TyrDC(EF_0634) 는 참조 세트에 포함되어 있으므로")
        print("     100% 에 가까워야 합니다. self-score 나 프로파일을 확인하십시오.")

    with open(LOG, "a") as lg:
        lg.write(f"\n[04_validate_controls] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  대조군 {len(todo)} · 통과 {ok} · 실패 {len(fail)}\n")
        lg.write(f"  self-score {self_score} · 임계 {CUTOFF}%\n")
        if v583:
            lg.write(f"  V583 자기대조 {v583:.2f}%\n")
        for n, e, g, p in fail:
            lg.write(f"  실패: {n} 기대={e} 실제={g} ({p:.2f}%)\n")

    if fail:
        print("\n" + "=" * 72)
        print("★ 여기서 멈춥니다. 전체 분석으로 넘어가지 마십시오.")
        print("=" * 72)
        print("  확인할 것:")
        print("   1) refs/tyrdc.hmm 과 refs/self_score.txt 가 같은 실행에서 나온 것인가")
        print("      (프로파일만 다시 만들고 self-score 를 갱신하지 않으면 척도가 어긋납니다)")
        print("   2) 그 유전체의 protein.faa 가 제대로 받아졌는가")
        print("   3) 실패한 것이 음성 대조인가 — 그렇다면 참조 세트가 너무 넓습니다")
        sys.exit(2)

    print("\n=== 4단계 통과 ===")
    print("다음:  bash 05_search_all.sh")


if __name__ == "__main__":
    main()
