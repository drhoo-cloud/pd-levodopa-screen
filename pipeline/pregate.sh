#!/usr/bin/env bash
set -uo pipefail
mkdir -p search/pregate logs
CARD=db/amr/protein_fasta_protein_homolog_model.fasta
VFDB=db/amr/vfdb.fas
echo "=== 프리게이트: Enterococcus 대조군 ==="
N=$(wc -l < panel/accessions_ent.txt); echo "  대상 $N 균주"

echo "  단백질 추출"
python - <<'PY'
keep=set(l.strip() for l in open('panel/accessions_ent.txt') if l.strip())
out=open('search/pregate/ent.faa','w'); w=False; n=0
for line in open('proteins/all_proteins.faa'):
    if line[0]=='>':
        w = line[1:].split('|')[0] in keep
        if w: n+=1
    if w: out.write(line)
out.close(); print(f'    {n} 서열')
PY

for DB in CARD VFDB; do
  [ "$DB" = CARD ] && F=$CARD || F=$VFDB
  echo "  $DB 검색"
  diamond makedb --in "$F" -d "search/pregate/${DB}" --quiet
  diamond blastp -q search/pregate/ent.faa -d "search/pregate/${DB}" \
    -o "search/pregate/${DB}_hits.tsv" --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --evalue 1e-10 --id 80 --query-cover 70 --max-target-seqs 5 --threads "$(nproc)" --quiet
  echo "    히트 $(wc -l < search/pregate/${DB}_hits.tsv) 건"
done

{ echo ""; echo "[pregate] $(date -Is)"; echo "  CARD $(wc -l < search/pregate/CARD_hits.tsv) · VFDB $(wc -l < search/pregate/VFDB_hits.tsv)"; echo "  기준 identity>=80 coverage>=70 evalue<=1e-10"; } | tee -a logs/run_log.txt
echo "=== 4단계 완료 ==="
