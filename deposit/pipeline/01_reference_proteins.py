#!/usr/bin/env python3
# ============================================================
# 01_reference_proteins.py  —  참조 단백질 세트 만들기
# ============================================================
# 실행:  python 01_reference_proteins.py
#
# 무엇을 하나
#   Gate 2 검색에 쓸 참조 단백질과 decoy 를 NCBI 에서 받아 FASTA 로 만듭니다.
#
# 왜 decoy 가 필요한가
#   TyrDC 는 group II 피리독살인산 탈탄산효소입니다.
#   같은 계열의 GadB(글루탐산), LdcA(라이신), OdcA(오르니틴) 과 서열이 상당히 닮았습니다.
#   TyrDC 만 놓고 검색하면 GadB 히트가 TyrDC 로 잘못 배정됩니다.
#   decoy 를 함께 넣고 "가장 높은 점수를 받은 계열" 로 배정해야 오배정을 막습니다.
#
# ★ 사람이 확인하는 단계가 들어 있습니다
#   받은 서열의 헤더를 눈으로 보고, 의도한 단백질이 맞는지 확인한 뒤 다음 단계로 갑니다.
# ============================================================

import os
import sys
import time
import datetime
import requests

OUT_DIR = "refs"
LOG = "logs/run_log.txt"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TODAY = datetime.date.today().isoformat()

# ------------------------------------------------------------
# 검색어 — accession 을 직접 적지 않고 질의로 받습니다.
#          accession 은 사람이 확인한 뒤 아래 VERIFIED 에 옮겨 적습니다.
# ------------------------------------------------------------
QUERIES = {
    # --- Gate 2 표적 ---
    "TARGET_tyrDC": [
        'tyrosine decarboxylase[Protein Name] AND "Enterococcus faecalis"[Organism]',
        'tyrosine decarboxylase[Protein Name] AND "Enterococcus faecium"[Organism]',
        'tyrosine decarboxylase[Protein Name] AND "Levilactobacillus brevis"[Organism]',
        'tyrosine decarboxylase[Protein Name] AND "Lactiplantibacillus plantarum"[Organism]',
    ],
    "TARGET_tyrP": [   # 타이로신 퍼미아제 — tyrDC 와 같은 오페론에 붙어 있습니다
        'tyrosine permease[Protein Name] AND "Enterococcus faecalis"[Organism]',
        'tyrosine permease[Protein Name] AND "Levilactobacillus brevis"[Organism]',
    ],
    "TARGET_aadc": [   # 방향족 아미노산 탈탄산효소 (넓은 계열)
        'aromatic amino acid decarboxylase[Protein Name] AND Lactobacillales[Organism]',
    ],
    # --- decoy: 같은 group II 계열이지만 표적이 아님 ---
    "DECOY_gadB": [
        'glutamate decarboxylase[Protein Name] AND Lactobacillales[Organism]',
    ],
    "DECOY_ldcA": [
        'lysine decarboxylase[Protein Name] AND Lactobacillales[Organism]',
    ],
    "DECOY_odcA": [
        'ornithine decarboxylase[Protein Name] AND Lactobacillales[Organism]',
    ],
    # --- 프리게이트: 히스티딘 탈탄산효소는 pyruvoyl 의존형이라 표적 계열이 아님 ---
    "PREGATE_hdcA": [
        'histidine decarboxylase[Protein Name] AND Lactobacillales[Organism]',
    ],
}

MAX_PER_QUERY = 5   # 질의당 최대 몇 건을 받을지


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
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    all_records = []          # (label, header, seq)
    summary = []

    print("=== 1단계: 참조 단백질 받기 ===")
    print("NCBI protein 데이터베이스에 질의합니다. 몇 분 걸립니다.\n")

    for label, terms in QUERIES.items():
        got = 0
        for term in terms:
            try:
                ids = esearch(term, MAX_PER_QUERY)
                time.sleep(0.34)          # NCBI 권장 간격 (초당 3회 이하)
                fasta = efetch(ids)
                time.sleep(0.34)
            except Exception as e:
                print(f"  ! {label}: {e}")
                continue

            # FASTA 를 레코드 단위로 쪼갭니다
            for block in fasta.split("\n>"):
                block = block.strip()
                if not block:
                    continue
                if not block.startswith(">"):
                    block = ">" + block
                lines = block.split("\n")
                header, seq = lines[0], "".join(lines[1:])
                if len(seq) < 100:        # 단편은 버립니다
                    continue
                all_records.append((label, header, seq))
                got += 1
        summary.append((label, got))
        print(f"  {label:16s} {got:3d} 건")

    # ------------------------------------------------------------
    # 파일로 쓰기 — 헤더 앞에 계열 라벨을 붙입니다.
    #   나중에 검색 결과에서 "이 히트가 어느 계열인지" 를 바로 알 수 있습니다.
    # ------------------------------------------------------------
    out_fa = os.path.join(OUT_DIR, "gate2_reference.faa")
    with open(out_fa, "w") as f:
        for label, header, seq in all_records:
            acc = header[1:].split()[0]
            desc = " ".join(header[1:].split()[1:])
            f.write(f">{label}|{acc} {desc}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    with open(LOG, "a") as lg:
        lg.write(f"\n[01_reference_proteins] {datetime.datetime.now().isoformat()}\n")
        for label, n in summary:
            lg.write(f"  {label}: {n}\n")
        lg.write(f"  총 {len(all_records)} 서열 → {out_fa}\n")

    print(f"\n총 {len(all_records)} 서열을 {out_fa} 에 저장했습니다.")
    print("\n" + "=" * 62)
    print("★ 여기서 멈추고 사람이 확인하십시오")
    print("=" * 62)
    print(f"1) 파일을 여십시오:   grep '>' {out_fa} | head -40")
    print("2) 확인할 것:")
    print("   · TARGET_tyrDC 에 정말 tyrosine decarboxylase 만 들어 있는가")
    print("   · DECOY_gadB 에 glutamate decarboxylase 가 들어 있는가")
    print("   · 'putative', 'hypothetical' 만 있는 항목은 지우십시오")
    print("   · 계열마다 최소 2건 이상 남아 있어야 합니다")
    print("3) 확인이 끝나면 다음 단계로:")
    print("   python 02_build_panel.py --taxa taxa_list.txt --out panel/")
    print()
    print("확인한 accession 목록은 원고 Supplementary Table S7 에 그대로 넣습니다.")


if __name__ == "__main__":
    main()
