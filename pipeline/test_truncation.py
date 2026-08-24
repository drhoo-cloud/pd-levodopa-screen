#!/usr/bin/env python3
# ============================================================
# test_truncation.py  —  변이 분류가 실제로 발동하는지 확인
# ============================================================
# 실행:  python test_truncation.py
#
# 무엇을 하나
#   참조 서열을 인위로 망가뜨린 합성 대조를 만들어
#   06_call_gates.py 의 분류가 그것들을 잡아내는지 봅니다.
#
#     TRUNC_xx   앞에서부터 xx% 만 남긴 절단 서열   → truncated 로 잡혀야 함
#     POINT_xx   무작위 위치 xx% 를 치환한 서열      → point-variant 로 잡혀야 함
#     INTACT     원본 그대로                         → intact 로 잡혀야 함
#
# ★ 왜 이 확인이 필요한가
#   어떤 분류가 실제 데이터에서 한 번도 안 나온다고 합시다.
#   그 계열이 정말 없는 것인지, 임계값 조합 때문에 원리상 나올 수 없는 것인지
#   구분할 방법이 없습니다. 합성 대조는 그것을 구분해 줍니다.
#
#   "truncated 가 데이터에 없다" 는 결론은
#   이 시험에서 truncated 가 나온 뒤에만 의미가 있습니다.
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점
# ------------------------------------------------------------
#   DIAMOND + identity/coverage/margin + decoy → hmmsearch + 정규화 점수
#   임계값은 06_call_gates.py 에서 직접 읽어옵니다.
#   두 파일에 같은 숫자를 두 번 적으면 반드시 어긋납니다.
# ============================================================

import os
import re
import sys
import random
import subprocess
import collections

random.seed(42)

REF = "refs/tyrdc_reference.faa"
HMM = "refs/tyrdc.hmm"
SELF_F = "refs/self_score.txt"
OUT = "search/trunc_test"
AA = "ACDEFGHIKLMNPQRSTVWY"


def load_thresholds():
    """06_call_gates.py 에서 임계값을 그대로 읽어옵니다."""
    src = open("06_call_gates.py").read()

    def grab(name, default):
        m = re.search(rf"^{name}\s*=\s*([0-9.]+)", src, re.M)
        return float(m.group(1)) if m else default

    return grab("CUTOFF", 50.0), grab("INTACT_MIN", 90.0), grab("TRUNC_RATIO", 0.80)


def read_fasta(path):
    seqs, lab, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if lab:
                seqs[lab] = "".join(buf)
            lab, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if lab:
        seqs[lab] = "".join(buf)
    return seqs


def main():
    for f in (REF, HMM, SELF_F):
        if not os.path.exists(f):
            sys.exit(f"파일이 없습니다: {f}  (01 단계를 먼저 실행하십시오)")

    CUTOFF, INTACT_MIN, TRUNC_RATIO = load_thresholds()
    self_score = float(open(SELF_F).read().strip())
    os.makedirs(OUT, exist_ok=True)

    seqs = read_fasta(REF)
    if not seqs:
        sys.exit("참조 서열이 비어 있습니다.")

    # 가장 긴 참조를 기준으로 삼습니다
    base = max(seqs, key=lambda k: len(seqs[k]))
    S = seqs[base]
    L = len(S)

    print("=== 변이 분류 발동 시험 ===")
    print(f"  기준 서열   {base}  ({L} aa)")
    print(f"  self-score  {self_score}")
    print(f"  임계값      present>{CUTOFF}%  intact>={INTACT_MIN}%  "
          f"절단<{TRUNC_RATIO}")
    print("  (임계값은 06_call_gates.py 에서 읽어왔습니다)\n")

    cases = []
    for frac in (0.95, 0.90, 0.85, 0.80, 0.70, 0.50):
        cases.append((f"TRUNC_{int(frac*100)}", S[:int(L * frac)]))
    for pct in (5, 10, 20, 40):
        s = list(S)
        for i in random.sample(range(L), int(L * pct / 100)):
            s[i] = random.choice(AA)
        cases.append((f"POINT_{pct}", "".join(s)))
    cases.append(("INTACT", S))

    q = os.path.join(OUT, "query.faa")
    with open(q, "w") as f:
        for n, s in cases:
            f.write(f">SYNTH_{n}|{n} synthetic control\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")
    print(f"  합성 대조 {len(cases)} 건 생성\n")

    dom = os.path.join(OUT, "hits.domtbl")
    subprocess.run(["hmmsearch", "--domtblout", dom, "-o", "/dev/null",
                    "-E", "1e-5", HMM, q], check=True)

    best = {}
    for line in open(dom):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 8:
            continue
        name, tlen, mlen, score = p[0], int(p[2]), int(p[5]), float(p[7])
        if name not in best or score > best[name]["score"]:
            best[name] = {"score": score, "tlen": tlen, "mlen": mlen}

    def classify(pct, tlen, mlen):
        if pct <= CUTOFF:
            return "absent", ""
        if mlen and tlen / mlen < TRUNC_RATIO:
            return "present", "truncated"
        if pct >= INTACT_MIN:
            return "present", "intact"
        return "present", "point-variant"

    print(f"  {'case':<12} {'score%':>8} {'len/model':>10}   판정 / 변이")
    print("  " + "-" * 58)
    seen = collections.Counter()
    for n, _ in cases:
        key = f"SYNTH_{n}|{n}"
        h = best.get(key)
        if not h:
            print(f"  {n:<12} {'히트 없음':>8}")
            seen["absent"] += 1
            continue
        pct = 100.0 * h["score"] / self_score
        ratio = h["tlen"] / max(h["mlen"], 1)
        call, var = classify(pct, h["tlen"], h["mlen"])
        print(f"  {n:<12} {pct:8.2f} {ratio:10.2f}   {call} / {var or '-'}")
        seen[var or call] += 1

    # ------------------------------------------------------------
    # 판정 — 세 분류가 모두 최소 한 번씩 나와야 합니다
    # ------------------------------------------------------------
    print("\n  분류별 발동 횟수")
    for k in ("intact", "point-variant", "truncated", "absent"):
        print(f"    {k:<15} {seen.get(k, 0)}")

    missing = [k for k in ("intact", "point-variant", "truncated")
               if seen.get(k, 0) == 0]
    print()
    if missing:
        print("  ★ 한 번도 발동하지 않은 분류: " + ", ".join(missing))
        print("    이 분류들은 현재 임계값 조합에서 원리상 나올 수 없습니다.")
        print("    실제 데이터에 없다는 결론을 내리기 전에 임계값을 다시 보십시오.")
        print("    특히 truncated 가 안 나오면 TRUNC_RATIO 가 너무 낮은 것입니다.")
    else:
        print("  세 분류가 모두 발동했습니다.")
        print("  실제 데이터에서 어떤 분류가 0건이라면, 그것은 진짜 0건입니다.")

    with open("logs/run_log.txt", "a") as lg:
        lg.write("\n[test_truncation] 합성 대조 %d건 · 발동 %s · 미발동 %s\n"
                 % (len(cases), dict(seen), missing or "없음"))


if __name__ == "__main__":
    main()
