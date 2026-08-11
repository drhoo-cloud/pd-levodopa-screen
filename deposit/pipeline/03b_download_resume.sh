#!/usr/bin/env bash
set -uo pipefail
ACC_FILE="${1:-panel/accessions_lp.txt}"
CHUNK=${CHUNK:-50}
MAX_RETRY=5
mkdir -p genomes proteins logs panel/chunks
TOTAL=$(wc -l < "$ACC_FILE")
echo "=== 3b 단계: 이어받기 ==="
echo "  전체 $TOTAL 건 · 묶음 $CHUNK"
find genomes -name "protein.faa" 2>/dev/null | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | sort -u > panel/_have.txt
HAVE=$(wc -l < panel/_have.txt); echo "  이미 확보 $HAVE 건"
sort -u "$ACC_FILE" > panel/_want.txt
comm -23 panel/_want.txt panel/_have.txt > panel/_todo.txt
TODO=$(wc -l < panel/_todo.txt); echo "  받을 것 $TODO 건"; echo
if [ "$TODO" -gt 0 ]; then
  rm -f panel/chunks/*
  split -l "$CHUNK" -d -a 3 panel/_todo.txt panel/chunks/c_
  for CH in panel/chunks/c_*; do
    N=$(wc -l < "$CH"); OK=0
    for TRY in 1 2 3 4 5; do
      printf "  %s (%d건) 시도 %d ... " "$(basename $CH)" "$N" "$TRY"
      if datasets download genome accession --inputfile "$CH" --include protein,gff3 --no-progressbar --filename "genomes/$(basename $CH).zip" >/dev/null 2>&1; then
        echo "성공"; OK=1; break
      fi
      W=$((TRY*15)); echo "실패 — ${W}초 후 재시도"
      rm -f "genomes/$(basename $CH).zip"; sleep $W
    done
    [ "$OK" -eq 0 ] && { echo "    건너뜀"; cat "$CH" >> logs/failed_accessions.txt; }
    sleep 2
  done
fi
echo "  압축 푸는 중"
cd genomes; for Z in *.zip; do [ -f "$Z" ] && unzip -oq "$Z" && rm -f "$Z"; done; cd ..
echo "  단백질 병합 중"
: > proteins/all_proteins.faa; : > logs/missing_protein_faa.txt
FOUND=0; MISSING=0
while read -r ACC; do
  P=$(find genomes -path "*${ACC}*" -name "protein.faa" 2>/dev/null | head -1)
  if [ -n "$P" ]; then
    awk -v acc="$ACC" '/^>/{print ">" acc "|" substr($0,2); next}{print}' "$P" >> proteins/all_proteins.faa
    FOUND=$((FOUND+1))
  else
    echo "$ACC" >> logs/missing_protein_faa.txt; MISSING=$((MISSING+1))
  fi
done < panel/_want.txt
NSEQ=$(grep -c '^>' proteins/all_proteins.faa 2>/dev/null || echo 0)
{ echo ""; echo "[03b] $(date -Is)"; echo "  요청 $TOTAL 확보 $FOUND 누락 $MISSING 서열 $NSEQ"; } | tee -a logs/run_log.txt
echo; echo "=== 3b 완료 ==="
echo "  확보 $FOUND / $TOTAL · 단백질 $NSEQ 서열"
echo "다음:  python 04_validate_controls.py"
