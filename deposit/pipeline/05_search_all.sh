#!/usr/bin/env bash
# ============================================================
# 05_search_all.sh  —  전체 패널 상동성 검색
# ============================================================
# 실행:  bash 05_search_all.sh
#
# 4단계(대조군)를 통과한 뒤에만 실행하십시오.
set -euo pipefail

[ -f search/controls.tsv ] || { echo "먼저 04_validate_controls.py 를 통과하십시오"; exit 1; }
echo "=== 5단계: 전체 검색 ==="

NSEQ=$(grep -c '^>' proteins/all_proteins.faa)
echo "  대상 단백질 $NSEQ 서열"

# Gate 2 — tyrDC/tdc + AADC + decoy
echo "  Gate 2 검색"
diamond blastp -q proteins/all_proteins.faa -d search/ref -o search/gate2_hits.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
  --evalue 1e-5 --max-target-seqs 25 --threads "$(nproc)" --quiet

G2=$(wc -l < search/gate2_hits.tsv)
echo "    히트 $G2 건"

# Gate 1 — CAZyme.  dbCAN3 이 있으면 그것을 쓰고, 없으면 안내만 합니다.
echo "  Gate 1 (CAZyme)"
if command -v run_dbcan >/dev/null 2>&1; then
  mkdir -p search/dbcan
  run_dbcan proteins/all_proteins.faa protein --out_dir search/dbcan --db_dir "$DBCAN_DB" \
    --tools hmmer --hmm_cpu "$(nproc)"
  echo "    dbCAN3 완료 → search/dbcan/overview.txt"
else
  cat <<'MSG'
    ! run_dbcan 이 없습니다. 아래로 설치하십시오.

        pip install dbcan
        mkdir -p db && cd db
        wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt
        hmmpress dbCAN-HMMdb-V12.txt
        cd ..
        export DBCAN_DB=$PWD/db

    설치 후 이 스크립트를 다시 실행하면 Gate 1 만 이어서 돌립니다.
MSG
fi

{
  echo ""
  echo "[05_search_all] $(date -Is)"
  echo "  단백질 $NSEQ · Gate2 히트 $G2"
  echo "  diamond: $(diamond --version 2>&1 | head -1)"
} | tee -a logs/run_log.txt

echo
echo "=== 5단계 완료 ==="
echo "다음:  python 06_call_gates.py"
