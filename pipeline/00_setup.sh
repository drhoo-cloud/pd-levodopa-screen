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
      ncbi-datasets-cli diamond blast hmmer mafft prodigal python=3.11
  pip install requests openpyxl pandas biopython statsmodels dbcan
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
  pip install --user requests openpyxl pandas biopython statsmodels dbcan

  # HMMER 와 MAFFT — 새 Gate 2 는 이 둘이 없으면 1단계에서 멈춥니다
  #   HMMER : 프로파일 생성(hmmbuild)과 검색(hmmsearch)
  #   MAFFT : 참조 서열 정렬. 프로파일의 입력입니다
  if ! command -v hmmbuild >/dev/null 2>&1 || ! command -v mafft >/dev/null 2>&1; then
    echo
    echo "  ! hmmbuild / mafft 가 없습니다. 단일 바이너리 배포본이 없어 자동 설치하지 않습니다."
    echo "    아래 중 하나로 설치하십시오."
    echo "        conda install -c bioconda hmmer mafft"
    echo "        sudo apt-get install -y hmmer mafft        # Debian/Ubuntu"
    echo "        sudo yum install -y hmmer mafft            # RHEL/CentOS"
  fi
fi

# ------------------------------------------------------------
# 2) 버전 기록 — 재현성의 출발점
# ------------------------------------------------------------
{
  echo "==============================================="
  echo "설치 일시: $(date -Is)"
  echo "-----------------------------------------------"
  for t in datasets diamond blastp hmmsearch hmmbuild mafft run_dbcan prodigal python3; do
    if command -v $t >/dev/null 2>&1; then
      printf "%-12s %s\n" "$t" "$($t --version 2>&1 | head -1)"
    else
      printf "%-12s (없음)\n" "$t"
    fi
  done
  echo "==============================================="
} | tee -a logs/run_log.txt

# ------------------------------------------------------------
# 3) 필수 도구 확인 — 없으면 여기서 멈춥니다
#    01 단계에 들어가서 실패하는 것보다 지금 아는 편이 낫습니다
# ------------------------------------------------------------
MISSING=""
for t in mafft hmmbuild hmmsearch diamond; do
  command -v "$t" >/dev/null 2>&1 || MISSING="$MISSING $t"
done
if [ -n "$MISSING" ]; then
  echo
  echo "★ 필수 도구가 없습니다:$MISSING"
  echo "  설치한 뒤 이 스크립트를 다시 실행하십시오."
  exit 1
fi

# ------------------------------------------------------------
# 4) dbCAN 데이터베이스 — 세 가지가 모두 있어야 합니다
#    HMMER 것만 받으면 --tools all 을 줘도 단일 method 로 돌아
#    GH family 가 과다 계상됩니다. 이번 개정의 원인이 된 지점입니다.
# ------------------------------------------------------------
if [ ! -d db ] || [ ! -f db/dbCAN-HMMdb-V12.txt.h3i ] \
   || [ ! -f db/CAZy.dmnd ] || [ ! -f db/dbCAN_sub.hmm.h3i ]; then
  cat <<'MSG'

! dbCAN 데이터베이스가 완전하지 않습니다. 세 가지가 모두 필요합니다.

    mkdir -p db && cd db
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt
    hmmpress dbCAN-HMMdb-V12.txt
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa
    diamond makedb --in CAZyDB.07262023.fa -d CAZy
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN_sub.hmm
    hmmpress dbCAN_sub.hmm
    cd ..
    export DBCAN_DB=$PWD/db

  HMMER 것만 받으면 Gate 1 이 단일 method 로 돌아갑니다.
MSG
else
  echo "dbCAN 데이터베이스 3종 확인 완료"
fi

echo
echo "=== 0단계 완료 ==="
echo "다음:  python 01_reference_proteins.py"
