#!/usr/bin/env python3
# ============================================================
# 01_reference_proteins.py  —  Gate 2 참조 세트와 프로파일 HMM 만들기
# ============================================================
# 실행:  python 01_reference_proteins.py
#
# 무엇을 하나
#   UniProtKB 에서 "검증된 세균 타이로신 탈탄산효소" 만 받아
#   프로파일 HMM 을 만들고, 정규화 분모가 되는 self-score 를 구합니다.
#
# ------------------------------------------------------------
# 왜 decoy 를 쓰지 않는가  (구 버전에서 바뀐 핵심)
# ------------------------------------------------------------
#   구 버전은 GadB(글루탐산), LdcA(라이신), OdcA(오르니틴) 을 decoy 로 함께 넣고
#   "가장 높은 점수를 받은 계열" 로 배정했습니다. 이 설계는 옳지 않습니다.
#
#   · 어떤 유전자가 있느냐 없느냐는, 비슷한 프로파일 여럿 중
#     무엇이 우연히 1등을 했는지로 정해지지 않습니다.
#   · 두 group II 프로파일 사이의 점수 차이(margin)는
#     표적 유전자가 있는지에 대해 아무것도 말해주지 않습니다.
#
#   새 설계는 질문을 하나로 좁힙니다.
#     "이 유전체에서 가장 좋은 단백질이, 검증된 세균 TyrDC 를 얼마나 닮았는가"
#   그리고 그것을 프로파일 자기점수 대비 백분율이라는 절대 척도로 답합니다.
#
#   ODC / HDC / GDC 는 levodopa 를 탈탄산하지 않으므로
#   참조 세트에서 '제외' 합니다. decoy 로 '포함' 하지 않습니다.
#
# ------------------------------------------------------------
# 산출물
# ------------------------------------------------------------
#   refs/tyrdc_reference.faa   참조 서열 (Zenodo 에 그대로 올립니다)
#   refs/tyrdc_reference.tsv   accession · 유전자명 · 생물종 · EC · 근거
#   refs/tyrdc.aln             MAFFT 정렬
#   refs/tyrdc.hmm             프로파일 HMM
#   refs/self_score.txt        정규화 분모 (이 값이 없으면 % 척도를 복원할 수 없습니다)
#   refs/TableS5_rows.tsv      Supplementary Table S5 에 붙여넣을 행
#   refs/pregate_hdc.faa       프리게이트용 히스티딘 탈탄산효소 참조
#
# ★ 사람이 확인하는 단계가 들어 있습니다
# ============================================================

import os
import sys
import csv
import json
import shutil
import hashlib
import subprocess
import datetime
import requests

OUT_DIR = "refs"
LOG = "logs/run_log.txt"
STREAM = "https://rest.uniprot.org/uniprotkb/stream"
TODAY = datetime.date.today().isoformat()

# ------------------------------------------------------------
# 질의 — accession 을 직접 적지 않습니다.
#         적어 두면 UniProt 이 개정될 때 조용히 낡습니다.
#         질의를 적어 두면 언제든 다시 만들 수 있습니다.
# ------------------------------------------------------------
#   gene / protein_name : locus 명으로 붙은 항목과 효소명으로 붙은 항목을 모두 잡습니다
#   ec:4.1.1.*          : 카르복시-리아제로 한정. 이름만 같고 활성이 다른 항목을 걸러냅니다
#   reviewed:true       : Swiss-Prot 만. 자동 주석 전이가 아니라 실험 근거가 있는 항목만
#   taxonomy_id:2       : 세균만. 식물·고세균 TyrDC 는 너무 멀어 프로파일을 흐립니다
# ------------------------------------------------------------
QUERY_TYRDC = (
    '((gene:tyrdc) OR (gene:tdc) OR (protein_name:"Tyrosine decarboxylase"))'
    ' AND (ec:4.1.1.*) AND (reviewed:true) AND (taxonomy_id:2)'
)

# 프리게이트용 — 히스티딘 탈탄산효소.
#   젖산균의 HDC 는 pyruvoyl 의존형이라 group II 계열이 아닙니다.
#   Gate 2 의 표적이 아니므로 여기서 따로 받아 프리게이트에서 씁니다.
QUERY_HDC = (
    '(protein_name:"Histidine decarboxylase")'
    ' AND (ec:4.1.1.22) AND (reviewed:true) AND (taxonomy_id:2)'
)

FIELDS = ("accession,id,protein_name,gene_names,organism_name,ec,"
          "length,sequence_version,protein_existence,annotation_score")

MIN_EXPECTED = 2      # 이보다 적게 나오면 질의가 잘못된 것으로 봅니다


# ------------------------------------------------------------
# UniProt 내려받기
# ------------------------------------------------------------
def fetch(query, fmt, fields=None):
    params = {"query": query, "format": fmt}
    if fields:
        params["fields"] = fields
    r = requests.get(STREAM, params=params, timeout=180)
    r.raise_for_status()
    return r.text, r.headers.get("x-uniprot-release", "UNKNOWN")


def sh(cmd):
    """외부 도구 실행. 실패하면 그 자리에서 멈춥니다."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"\n실행 실패: {' '.join(cmd)}\n{p.stderr[:800]}")
    return p.stdout


def need(tool, howto):
    if not shutil.which(tool):
        sys.exit(f"\n{tool} 이 없습니다.\n  설치:  {howto}")


# ------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    need("mafft", "conda install -c bioconda mafft")
    need("hmmbuild", "conda install -c bioconda hmmer")

    print("=== 1단계: Gate 2 참조 세트와 프로파일 만들기 ===\n")

    # --- 1) 표적 서열 ---------------------------------------
    print("  UniProtKB 질의")
    print(f"    {QUERY_TYRDC}")
    fasta, release = fetch(QUERY_TYRDC, "fasta")
    tsv, _ = fetch(QUERY_TYRDC, "tsv", FIELDS)

    faa = os.path.join(OUT_DIR, "tyrdc_reference.faa")
    open(faa, "w").write(fasta)
    open(os.path.join(OUT_DIR, "tyrdc_reference.tsv"), "w").write(tsv)

    rows = list(csv.DictReader(tsv.splitlines(), delimiter="\t"))
    n = len(rows)
    print(f"    UniProt 릴리스 {release} · {n} 건\n")

    if n < MIN_EXPECTED:
        sys.exit(f"★ {n} 건뿐입니다. 질의가 의도대로 걸리지 않았습니다.\n"
                 f"  UniProt 웹에서 같은 질의를 먼저 눈으로 확인하십시오.")

    # 받은 것 보여주기 — 여기서 사람이 한 번 봅니다
    def first(r, *keys):
        for k in keys:
            if r.get(k):
                return r[k]
        return ""

    print("  받은 서열")
    print("  " + "-" * 72)
    table = []
    for r in rows:
        acc = first(r, "Entry", "accession")
        gene = (first(r, "Gene Names", "gene_names").split() or [""])[0]
        org = first(r, "Organism", "organism_name").split("(")[0].strip()
        ec = first(r, "EC number", "ec")
        pe = first(r, "Protein existence", "protein_existence")
        score = first(r, "Annotation score", "annotation_score")
        length = first(r, "Length", "length")
        table.append([acc, gene, " ".join(org.split()[:2]), ec, length, pe, score])
        print(f"  {acc:12s} {gene:10s} {org[:34]:34s} {ec:9s} {length:>4s} aa")
    print("  " + "-" * 72 + "\n")

    # Table S5 에 붙여넣을 행
    with open(os.path.join(OUT_DIR, "TableS5_rows.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Accession", "Gene", "Organism", "EC",
                    "Length", "Protein existence", "Annotation score"])
        w.writerows(sorted(table, key=lambda x: x[2]))

    # --- 2) 프리게이트용 HDC --------------------------------
    print("  프리게이트용 히스티딘 탈탄산효소 받기")
    try:
        hdc, _ = fetch(QUERY_HDC, "fasta")
        open(os.path.join(OUT_DIR, "pregate_hdc.faa"), "w").write(hdc)
        print(f"    {hdc.count('>')} 건 → refs/pregate_hdc.faa\n")
    except Exception as e:
        print(f"    ! 실패 (프리게이트는 나중에 따로 돌릴 수 있습니다): {e}\n")

    # --- 3) 정렬과 프로파일 --------------------------------
    aln = os.path.join(OUT_DIR, "tyrdc.aln")
    hmm = os.path.join(OUT_DIR, "tyrdc.hmm")

    print("  MAFFT 정렬")
    open(aln, "w").write(sh(["mafft", "--auto", faa]))

    print("  hmmbuild 프로파일 생성")
    build_log = sh(["hmmbuild", hmm, aln])

    mlen = ""
    for line in build_log.splitlines():
        p = line.split()
        if len(p) >= 6 and p[0].isdigit():
            mlen = p[3]

    # --- 4) self-score ------------------------------------
    #   프로파일을 '자기가 만들어진 서열' 에 되돌려 검색합니다.
    #   그때 나오는 최고 점수가 정규화의 분모입니다.
    #   이 값이 없으면 논문에 적힌 93.0% 같은 수치를 아무도 복원할 수 없습니다.
    print("  self-score 산출")
    tbl = os.path.join(OUT_DIR, "self.tbl")
    sh(["hmmsearch", "--tblout", tbl, "-o", os.path.join(OUT_DIR, "self.txt"),
        hmm, faa])

    scores = []
    for line in open(tbl):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) > 5:
            scores.append((p[0], float(p[5])))
    if not scores:
        sys.exit("★ 자기 히트가 하나도 없습니다. 프로파일이 제대로 만들어지지 않았습니다.")

    self_score = max(s for _, s in scores)
    open(os.path.join(OUT_DIR, "self_score.txt"), "w").write(f"{self_score}\n")
    md5 = hashlib.md5(open(hmm, "rb").read()).hexdigest()

    print("  " + "-" * 72)
    for name, s in sorted(scores, key=lambda x: -x[1]):
        print(f"  {name:32s} {s:8.1f}  {100*s/self_score:6.2f}%")
    print("  " + "-" * 72)
    print(f"  self-score = {self_score}   (정규화 분모)")
    print(f"  md5        = {md5}\n")

    # --- 5) 기록 -------------------------------------------
    meta = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "uniprot_release": release,
        "retrieval_date": TODAY,
        "query_tyrdc": QUERY_TYRDC,
        "query_hdc": QUERY_HDC,
        "n_entries": n,
        "accessions": [t[0] for t in table],
        "aligner": sh(["mafft", "--version"]).strip() or "mafft",
        "model_length": mlen,
        "profile_md5": md5,
        "self_score": self_score,
        "cutoff_percent": 50.0,
    }
    json.dump(meta, open(os.path.join(OUT_DIR, "gate2_profile.json"), "w"),
              indent=2, ensure_ascii=False)

    with open(LOG, "a") as lg:
        lg.write(f"\n[01_reference_proteins] {meta['generated']}\n")
        lg.write(f"  UniProt 릴리스 {release} · 검색일 {TODAY}\n")
        lg.write(f"  질의 {QUERY_TYRDC}\n")
        lg.write(f"  서열 {n} 건: {', '.join(meta['accessions'])}\n")
        lg.write(f"  모델 길이 {mlen} · md5 {md5}\n")
        lg.write(f"  self-score {self_score} · 임계 50%\n")

    # --- 6) 사람이 확인 ------------------------------------
    print("=" * 74)
    print("★ 여기서 멈추고 확인하십시오")
    print("=" * 74)
    print("1) 위 표에 tyrosine decarboxylase 만 있는지 보십시오.")
    print("   permease, aromatic amino acid decarboxylase 가 섞였다면 질의를 좁히십시오.")
    print("2) Protein existence 가 모두 'protein level' 인지 보십시오.")
    print("3) self-score 를 논문 본문에 보고한 값과 대조하십시오.")
    print("   참조 서열의 정규화 점수가 종별 보고값을 재현해야 합니다.")
    print("4) refs/TableS5_rows.tsv 를 Supplementary Table S5 에 넣으십시오.")
    print("   (구 버전 주석은 Table S7 이라고 되어 있었습니다. S5 가 맞습니다.)")
    print()
    print("다음:  python 02_build_panel.py --taxa taxa_list.txt --out panel/")


if __name__ == "__main__":
    main()
