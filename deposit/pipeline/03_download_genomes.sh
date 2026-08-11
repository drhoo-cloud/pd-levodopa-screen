#!/usr/bin/env bash
# ============================================================
# 03_download_genomes.sh  —  유전체와 단백질 서열 내려받기
# ============================================================
# 실행:  bash 03_download_genomes.sh panel/accessions_all.txt
set -euo pipefail

ACC_FILE="${1:-panel/accessions_all.txt}"
[ -f "$ACC_FILE" ] || { echo "파일이 없습니다: $ACC_FILE"; exit 1; }
N=$(wc -l < "$ACC_FILE")
echo "=== 3단계: 유전체 $N 건 내려받기 ==="

mkdir -p genomes proteins logs

# 한 번에 너무 많이 요청하면 끊깁니다. 200건씩 나눠 받습니다.
split -l 200 -d "$ACC_FILE" panel/_chunk_
for CH in panel/_chunk_*; do
  echo "  → $CH ($(wc -l < "$CH")건)"
  datasets download genome accession --inputfile "$CH" \
      --include genome,protein,gff3 \
      --filename "genomes/$(basename "$CH").zip" || {
        echo "    ! 실패 — 재시도"; sleep 10
        datasets download genome accession --inputfile "$CH" \
            --include genome,protein,gff3 \
            --filename "genomes/$(basename "$CH").zip"; }
  sleep 2
done

echo
echo "압축 풀기"
cd genomes
for Z in *.zip; do unzip -oq "$Z" && rm -f "$Z"; done
cd ..

# 단백질 서열을 한 파일로 모읍니다. 헤더 앞에 accession 을 붙여 출처를 잃지 않습니다.
echo "단백질 서열 병합"
: > proteins/all_proteins.faa
FOUND=0; MISSING=0
while read -r ACC; do
  P=$(find genomes -path "*${ACC}*" -name "protein.faa" | head -1)
  if [ -n "$P" ]; then
    awk -v acc="$ACC" '/^>/{print ">" acc "|" substr($0,2); next}{print}' "$P" \
      >> proteins/all_proteins.faa
    FOUND=$((FOUND+1))
  else
    echo "$ACC" >> logs/missing_protein_faa.txt
    MISSING=$((MISSING+1))
  fi
done < "$ACC_FILE"

TOTAL_SEQ=$(grep -c '^>' proteins/all_proteins.faa || echo 0)
{
  echo ""
  echo "[03_download_genomes] $(date -Is)"
  echo "  요청 $N · 단백질 확보 $FOUND · 누락 $MISSING"
  echo "  총 단백질 서열 $TOTAL_SEQ"
} | tee -a logs/run_log.txt

echo
echo "=== 3단계 완료 ==="
echo "  단백질 $TOTAL_SEQ 서열 → proteins/all_proteins.faa"
[ "$MISSING" -gt 0 ] && echo "  ★ 누락 $MISSING 건은 logs/missing_protein_faa.txt 를 보십시오"
echo "다음:  python 04_validate_controls.py"
