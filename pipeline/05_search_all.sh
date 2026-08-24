#!/usr/bin/env bash
# ============================================================
# 05_search_all.sh  —  전체 패널 검색
# ============================================================
# 실행:  bash 05_search_all.sh
#
# 4단계(대조군)를 통과한 뒤에만 실행하십시오.
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점
# ------------------------------------------------------------
#   Gate 2 : DIAMOND blastp (참조 서열 + decoy) → hmmsearch (프로파일 HMM)
#            identity 와 계열 간 margin 으로 판정하던 것을,
#            프로파일 자기점수 대비 백분율 하나로 바꿉니다.
#
#   Gate 1 : run_dbcan --tools hmmer → gate1_cazyme.sh (세 방법 중 둘 이상 일치)
#            여기서 직접 부르지 않고 별도 스크립트로 분리했습니다.
# ------------------------------------------------------------
set -euo pipefail

PROT=proteins/all_proteins.faa
HMM=refs/tyrdc.hmm
EVALUE=1e-5

[ -f search/controls.tsv ] || { echo "먼저 04_validate_controls.py 를 통과하십시오"; exit 1; }
[ -f "$HMM" ] || { echo "프로파일이 없습니다: $HMM  (01 단계를 먼저 실행하십시오)"; exit 1; }
[ -f refs/self_score.txt ] || { echo "self-score 가 없습니다: refs/self_score.txt"; exit 1; }

mkdir -p search logs
SELF=$(tr -d '[:space:]' < refs/self_score.txt)
NSEQ=$(grep -c '^>' "$PROT")

echo "=== 5단계: 전체 검색 ==="
echo "  대상 단백질 $NSEQ 서열"
echo "  프로파일    $HMM"
echo "  self-score  $SELF   (정규화 분모)"
echo "  E-value     $EVALUE"
echo "  시작        $(date -Is)"

# ------------------------------------------------------------
# Gate 2 — 프로파일 검색
#   --domtblout 을 씁니다. tblout 에는 표적 단백질 길이가 없어서
#   절단 여부(truncated)를 판정할 수 없습니다.
# ------------------------------------------------------------
hmmsearch --cpu "$(nproc)" -E "$EVALUE" \
  --tblout    search/gate2_hits.tbl \
  --domtblout search/gate2_hits.domtbl \
  -o          search/gate2_hits.txt \
  "$HMM" "$PROT"

G2=$(grep -vc '^#' search/gate2_hits.tbl || true)
echo "  검색 완료  $(date -Is)"
echo "  히트 $G2 건"

# ------------------------------------------------------------
# Gate 1 — 별도 스크립트
# ------------------------------------------------------------
echo
echo "  Gate 1 (CAZyme) 은 별도 스크립트로 돌립니다:"
echo "      bash gate1_cazyme.sh"
echo "  구 버전은 이 자리에서 run_dbcan --tools hmmer 를 불렀습니다."
echo "  단일 method 는 GH family 를 과다 계상하므로 더 이상 쓰지 않습니다."

{
  echo ""
  echo "[05_search_all] $(date -Is)"
  echo "  단백질 $NSEQ · Gate2 히트 $G2"
  echo "  hmmsearch: $(hmmsearch -h | head -2 | tail -1)"
  echo "  프로파일 $HMM · self-score $SELF · E-value $EVALUE"
} | tee -a logs/run_log.txt

echo
echo "=== 5단계 완료 ==="
echo "다음:  bash gate1_cazyme.sh  →  python 06_call_gates.py"
