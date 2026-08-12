#!/usr/bin/env python3
# ============================================================
# 01b_topup_references.py  —  건수가 적은 계열만 보충
# ============================================================
# 실행:  python 01b_topup_references.py
#
# 언제 쓰나
#   01 단계 결과에서 어떤 계열이 2건 미만일 때 씁니다.
#   특히 DECOY 계열이 얇으면 표적이 아닌 효소가 표적으로 잘못 배정됩니다.
#
# 무엇을 하나
#   검색어를 넓혀 부족한 계열만 더 받아 기존 파일에 덧붙입니다.
#   이미 있는 accession 은 건너뜁니다.
# ============================================================

import os
import re
import time
import datetime
import requests
from collections import Counter

FAA = "refs/gate2_reference.faa"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MIN_PER_LABEL = 5          # 계열마다 이 수를 목표로 채웁니다

# 01 단계보다 넓은 검색어. [Protein Name] 제한을 풀고 동의어를 넣었습니다.
WIDE = {
    "DECOY_ldcA": [
        'lysine decarboxylase AND Lactobacillales[Organism]',
        'lysine decarboxylase AND Firmicutes[Organism] AND 400:800[SLEN]',
        'ldcA AND Bacteria[Organism] AND 400:800[SLEN]',
        '"lysine/ornithine decarboxylase" AND Bacteria[Organism]',
    ],
    "DECOY_odcA": [
        'ornithine decarboxylase AND Lactobacillales[Organism]',
        'ornithine decarboxylase AND Firmicutes[Organism] AND 400:800[SLEN]',
    ],
    "DECOY_gadB": [
        'glutamate decarboxylase AND Lactobacillales[Organism]',
        'gadB AND Lactobacillales[Organism]',
    ],
    "TARGET_tyrP": [
        'tyrosine permease AND Lactobacillales[Organism]',
        'tyrP AND Lactobacillales[Organism]',
        'tyrosine permease AND Enterococcus[Organism]',
        'amino acid permease AND Enterococcus faecalis[Organism] AND 400:600[SLEN]',
    ],
    "TARGET_aadc": [
        'aromatic-L-amino-acid decarboxylase AND Bacteria[Organism]',
        'aromatic amino acid decarboxylase AND Firmicutes[Organism]',
        'group II decarboxylase AND Lactobacillales[Organism]',
    ],
    "TARGET_tyrDC": [
        'tyrosine decarboxylase AND Lactobacillales[Organism]',
        'tdcA AND Bacteria[Organism] AND 500:700[SLEN]',
    ],
    "PREGATE_hdcA": [
        'histidine decarboxylase AND Lactobacillales[Organism]',
        'hdcA AND Bacteria[Organism]',
    ],
}


def esearch(term, retmax):
    r = requests.get(f"{EUTILS}/esearch.fcgi", timeout=60, params={
        "db": "protein", "term": term, "retmax": retmax, "retmode": "json"})
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def efetch(ids):
    if not ids:
        return ""
    r = requests.get(f"{EUTILS}/efetch.fcgi", timeout=120, params={
        "db": "protein", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"})
    r.raise_for_status()
    return r.text


def main():
    if not os.path.exists(FAA):
        raise SystemExit(f"파일이 없습니다: {FAA}  (01 단계를 먼저 실행하십시오)")

    # 현재 상태 확인
    have = Counter()
    seen_acc = set()
    for line in open(FAA):
        if line.startswith(">"):
            head = line[1:].strip()
            label = head.split("|")[0]
            acc = head.split("|")[1].split()[0] if "|" in head else ""
            have[label] += 1
            seen_acc.add(acc)

    print("=== 1b 단계: 부족한 계열 보충 ===")
    print("현재 상태")
    for lab in sorted(WIDE):
        mark = "  " if have[lab] >= 2 else "★ "
        print(f"  {mark}{lab:16s} {have[lab]:3d} 건")

    thin = [l for l in WIDE if have[l] < MIN_PER_LABEL]
    if not thin:
        print("\n모든 계열이 충분합니다. 다음 단계로 가십시오.")
        return

    print(f"\n보충할 계열: {', '.join(thin)}\n")

    added = Counter()
    with open(FAA, "a") as out:
        for lab in thin:
            need = MIN_PER_LABEL - have[lab]
            for term in WIDE[lab]:
                if added[lab] >= need:
                    break
                try:
                    ids = esearch(term, 10)
                    time.sleep(0.34)
                    fasta = efetch(ids)
                    time.sleep(0.34)
                except Exception as e:
                    print(f"  ! {lab}: {e}")
                    continue

                for block in fasta.split("\n>"):
                    if added[lab] >= need:
                        break
                    block = block.strip()
                    if not block:
                        continue
                    if not block.startswith(">"):
                        block = ">" + block
                    lines = block.split("\n")
                    header, seq = lines[0], "".join(lines[1:])
                    acc = header[1:].split()[0]
                    desc = " ".join(header[1:].split()[1:])
                    if acc in seen_acc or len(seq) < 200:
                        continue
                    # 설명에 hypothetical / unnamed 만 있으면 버립니다
                    if re.search(r"hypothetical|unnamed|uncharacteri", desc, re.I):
                        continue
                    seen_acc.add(acc)
                    out.write(f">{lab}|{acc} {desc}\n")
                    for i in range(0, len(seq), 60):
                        out.write(seq[i:i+60] + "\n")
                    added[lab] += 1

    print("보충 결과")
    for lab in thin:
        print(f"  {lab:16s} +{added[lab]} → {have[lab] + added[lab]} 건")

    with open("logs/run_log.txt", "a") as lg:
        lg.write(f"\n[01b_topup_references] {datetime.datetime.now().isoformat()}\n")
        for lab in thin:
            lg.write(f"  {lab}: {have[lab]} -> {have[lab] + added[lab]}\n")

    print("\n다시 확인하십시오")
    print("  grep '>' refs/gate2_reference.faa | cut -d'|' -f1 | sort | uniq -c")
    print("  grep 'DECOY_ldcA' refs/gate2_reference.faa")


if __name__ == "__main__":
    main()
