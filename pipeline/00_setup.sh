#!/usr/bin/env bash
# ============================================================
# 00_setup.sh  —  도구 설치와 폴더 만들기 (한 번만 실행)
# ============================================================
# 실행:  bash 00_setup.sh
set -euo pipefail

echo "=== 0단계: 환경 준비 ==="
mkdir -p refs panel genomes proteins search results logs
echo "폴더 생성 완료: refs panel genomes proteins search results logs"
echo

# ------------------------------------------------------------
# 1) conda 가 있으면 conda 로, 없으면 pip + 바이너리로
# ------------------------------------------------------------
if command -v conda >/dev/null 2>&1; then
  echo "conda 를 찾았습니다. conda 로 설치합니다."
  conda install -y -c conda-forge -c bioconda \
      ncbi-datasets-cli diamond blast hmmer prodigal python=3.11
  pip install requests openpyxl pandas biopython statsmodels
else
  echo "conda 가 없습니다. 단일 바이너리로 설치합니다."
  mkdir -p bin && cd bin

  # NCBI datasets — 유전체 목록 조회와 다운로드에 씁니다
  if [ ! -f datasets ]; then
    curl -sSL -o datasets \
      'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
    curl -sSL -o dataformat \
      'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
    chmod +x datasets dataformat
  fi

  # DIAMOND — 단백질 상동성 검색. BLAST 보다 수백 배 빠릅니다
  if [ ! -f diamond ]; then
    curl -sSL https://github.com/bbuchfink/diamond/releases/latest/download/diamond-linux64.tar.gz \
      | tar xz diamond
    chmod +x diamond
  fi

  cd ..
  export PATH="$PWD/bin:$PATH"
  echo "export PATH=\"$PWD/bin:\$PATH\"" >> ~/.bashrc
  pip install --user requests openpyxl pandas biopython statsmodels
fi

# ------------------------------------------------------------
# 2) 버전 기록 — 재현성의 출발점
# ------------------------------------------------------------
{
  echo "==============================================="
  echo "설치 일시: $(date -Is)"
  echo "-----------------------------------------------"
  for t in datasets diamond blastp hmmsearch prodigal python3; do
    if command -v $t >/dev/null 2>&1; then
      printf "%-12s %s\n" "$t" "$($t --version 2>&1 | head -1)"
    else
      printf "%-12s (없음)\n" "$t"
    fi
  done
  echo "==============================================="
} | tee -a logs/run_log.txt

echo
echo "=== 0단계 완료 ==="
echo "다음:  python 01_reference_proteins.py"
