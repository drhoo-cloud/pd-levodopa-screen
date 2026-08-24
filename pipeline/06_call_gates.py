#!/usr/bin/env python3
# ============================================================
# 06_call_gates.py  —  유전체별 판정
# ============================================================
# 실행:  python 06_call_gates.py
#
# 무엇을 하나
#   Gate 2 : 유전체마다 정규화 프로파일 점수를 구하고 present / absent 를 정합니다.
#   Gate 1 : GH family 개수를 붙입니다. 통과/탈락을 매기지 않습니다.
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점  ★ 세 가지입니다
# ------------------------------------------------------------
#  (1) decoy 최고점 배정을 없앴습니다.
#      구 버전은 최상위 계열이 decoy 면 absent 로 처리했습니다.
#      어떤 유전자가 있느냐는 비슷한 프로파일 여럿의 순위 경쟁으로 정할 문제가 아닙니다.
#      이제 질문은 하나입니다 — 검증된 TyrDC 를 얼마나 닮았는가.
#
#  (2) ambiguous 범주를 없앴습니다.
#      구 버전은 identity/coverage 중간대와 margin 부족을 ambiguous 로 뒀습니다.
#      정규화 점수로 바꾸니 분포가 완전히 둘로 갈립니다.
#      5,373 유전체에서 음성은 0–4.3%, 양성은 61.5–100% 이고
#      그 사이에는 단 한 건(19.5%)뿐입니다. 중간 범주가 생기지 않습니다.
#
#      다만 임계 근처는 계속 표시합니다. NEAR_LO–NEAR_HI 구간에 들어오면
#      near_cutoff 열에 1 을 적어 효소 assay 로 회부하도록 합니다.
#      판정을 흐리지 않으면서 경고는 남기는 방식입니다.
#
#  (3) GH_MIN = 15 를 삭제했습니다.
#      원고는 Gate 1 에 대해 "No pass threshold is set" 을 세 곳에서 말합니다.
#      코드에 15 라는 기준이 남아 있으면 논문과 저장소가 어긋납니다.
#      Gate 1 의 출력은 판정이 아니라 분포상의 위치입니다.
# ============================================================

import os
import csv
import datetime
from collections import defaultdict

DOMTBL = "search/gate2_hits.domtbl"
SELF_F = "refs/self_score.txt"
DBCAN = "search/dbcan/overview_consensus.txt"
DBCAN_OLD = "search/dbcan/overview.txt"
PROT = "proteins/all_proteins.faa"
OUT = "results/per_genome_calls.tsv"
LOG = "logs/run_log.txt"

# ------------------------------------------------------------
# 사전 확정 임계값 — 04 단계와 같은 값이어야 합니다
# ------------------------------------------------------------
CUTOFF = 50.0            # 정규화 점수 이 값 초과면 present
NEAR_LO, NEAR_HI = 40.0, 60.0   # 임계 근처 — 판정은 하되 assay 회부 표시
INTACT_MIN = 90.0        # 이 이상이면 intact 로 봅니다
TRUNC_RATIO = 0.80       # 모델 길이의 80% 미만이면 절단


def load_self_score(path):
    if not os.path.exists(path):
        raise SystemExit(f"파일이 없습니다: {path}  (01 단계를 먼저 실행하십시오)\n"
                         f"  self-score 가 없으면 % 척도를 만들 수 없습니다.")
    v = float(open(path).read().strip())
    if v <= 0:
        raise SystemExit(f"self-score 가 이상합니다: {v}")
    return v


def load_hits(path):
    """hmmsearch --domtblout 을 단백질 단위로 읽습니다.

    열 위치 (1-based)
        1  target name   검색 대상 단백질
        3  tlen          그 단백질의 길이
        6  qlen          프로파일(모델) 길이
        8  score         full-sequence bit score   ← 이 값을 씁니다
    """
    best = {}
    if not os.path.exists(path):
        raise SystemExit(f"파일이 없습니다: {path}  (05 단계를 먼저 실행하십시오)")
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 8:
            continue
        name, tlen, mlen, score = p[0], int(p[2]), int(p[5]), float(p[7])
        cur = best.get(name)
        if cur is None or score > cur["score"]:
            best[name] = {"score": score, "tlen": tlen, "mlen": mlen}
    return best


def classify(pct, tlen, mlen):
    """정규화 점수 하나로 판정합니다. 계열 비교도 margin 도 쓰지 않습니다."""
    if pct <= CUTOFF:
        return "absent", ""
    if mlen and tlen / mlen < TRUNC_RATIO:
        # 유전자는 있으나 뚜렷하게 짧습니다 — 활성 확인이 필요합니다
        return "present", "truncated"
    if pct >= INTACT_MIN:
        return "present", "intact"
    # 온전한 길이지만 참조에서 상당히 떨어져 있습니다.
    # levodopa 를 실제로 탈탄산하는지는 서열만으로 말할 수 없습니다.
    return "present", "point-variant"


def load_gh_counts():
    """GH family 개수. gate1_cazyme.sh 의 consensus 출력을 우선합니다."""
    if os.path.exists(DBCAN):
        counts = {}
        with open(DBCAN) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                counts[row["accession"].split("|")[0]] = int(row["n_GH_families"])
        return counts, "consensus"

    if os.path.exists(DBCAN_OLD):
        print("  ! consensus 출력이 없어 구 overview.txt 를 씁니다.")
        print("    이 값은 단일 method 결과라 GH family 가 과다 계상됩니다.")
        print("    bash gate1_cazyme.sh 를 먼저 돌리십시오.")
        counts = defaultdict(set)
        with open(DBCAN_OLD) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                acc = p[0].split("|")[0]
                for field in p[1:]:
                    for tok in str(field).replace("+", " ").replace(";", " ").split():
                        if tok.startswith("GH"):
                            counts[acc].add(tok.split("_")[0].split("(")[0])
        return {k: len(v) for k, v in counts.items()}, "single-method(구버전)"

    return {}, "미실행"


def main():
    os.makedirs("results", exist_ok=True)

    self_score = load_self_score(SELF_F)
    best = load_hits(DOMTBL)
    gh, gh_source = load_gh_counts()

    # 히트가 없는 유전체도 분모에 포함합니다 (명백한 음성)
    all_acc = set()
    with open(PROT) as f:
        for line in f:
            if line.startswith(">"):
                all_acc.add(line[1:].split("|")[0])

    print("=== 6단계: 판정 ===")
    print(f"  전체 유전체 {len(all_acc)} 건 (히트 유무 무관)")
    print(f"  self-score  {self_score}  ·  임계 {CUTOFF}%")
    print(f"  Gate 1 출처 {gh_source}")

    # 단백질 단위 판정 → 유전체 단위로 모읍니다
    genome = defaultdict(lambda: {"hits": []})
    for name, h in best.items():
        acc = name.split("|")[0]
        pct = 100.0 * h["score"] / self_score
        call, variant = classify(pct, h["tlen"], h["mlen"])
        genome[acc]["hits"].append((name, pct, call, variant))

    rows, tally = [], defaultdict(int)
    near = 0
    for acc in sorted(all_acc):
        hits = genome[acc]["hits"]
        top_pct = max((p for _, p, _, _ in hits), default=0.0)
        present = [(n, p, v) for n, p, c, v in hits if c == "present"]

        if present:
            g2 = "present"
            variants = ";".join(sorted({v for _, _, v in present}))
            loci = ";".join(n for n, _, _ in present)
        else:
            g2, variants, loci = "absent", "", ""

        is_near = 1 if NEAR_LO <= top_pct <= NEAR_HI else 0
        near += is_near
        tally[g2] += 1

        rows.append({
            "assembly_accession": acc,
            "gate2_call": g2,
            "gate2_score_pct": f"{top_pct:.2f}",
            "gate2_variant": variants,
            "gate2_near_cutoff": is_near,
            "gate2_loci": loci,
            "gate1_gh_families": gh.get(acc, ""),
        })

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\n  Gate 2 present  {tally['present']:5d}")
    print(f"  Gate 2 absent   {tally['absent']:5d}")
    print(f"  임계 근처       {near:5d}   ({NEAR_LO}–{NEAR_HI}%, assay 회부 대상)")

    # ------------------------------------------------------------
    # 분포가 정말 둘로 갈리는지 확인합니다.
    # 중간대에 유전체가 몰려 있다면 임계값 위치가 결과를 좌우한다는 뜻이고,
    # 그때는 프로파일이나 참조 세트를 다시 봐야 합니다.
    # ------------------------------------------------------------
    pcts = [float(r["gate2_score_pct"]) for r in rows]
    below = [p for p in pcts if p <= CUTOFF]
    above = [p for p in pcts if p > CUTOFF]
    print("\n  점수 분포")
    if below:
        print(f"    음성쪽 최대 {max(below):6.2f}%")
    if above:
        print(f"    양성쪽 최소 {min(above):6.2f}%")
    if below and above:
        gap = min(above) - max(below)
        print(f"    간격        {gap:6.2f}%p", end="")
        print("   → 분포가 이분화되어 있습니다" if gap > 20
              else "   ★ 간격이 좁습니다. 임계값이 결과를 좌우할 수 있습니다")

    if gh:
        counts = sorted(v for v in gh.values())
        mid = counts[len(counts) // 2]
        print(f"\n  Gate 1 GH family 중앙값 {mid} · 범위 {min(counts)}–{max(counts)}")
        print("    통과/탈락은 매기지 않습니다. 분포상의 위치가 출력입니다.")
    else:
        print("\n  Gate 1 미실행 — bash gate1_cazyme.sh 를 돌리십시오")

    with open(LOG, "a") as lg:
        lg.write(f"\n[06_call_gates] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  유전체 {len(rows)} · present {tally['present']} "
                 f"· absent {tally['absent']} · 임계근처 {near}\n")
        lg.write(f"  self-score {self_score} · 임계 {CUTOFF}% "
                 f"· intact>={INTACT_MIN}% · 절단<{TRUNC_RATIO}\n")
        lg.write(f"  Gate 1 출처 {gh_source} · pass threshold 없음\n")

    print(f"\n  → {OUT}")
    print("다음:  python 07_summarize.py")


if __name__ == "__main__":
    main()
