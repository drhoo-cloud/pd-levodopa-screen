#!/usr/bin/env bash
set -uo pipefail
DB=db/dbCAN.hmm
[ -f "$DB" ] || { echo "DB 없음: $DB"; exit 1; }
[ -f "$DB.h3i" ] || { echo "hmmpress 필요: cd db && hmmpress dbCAN.hmm && cd .."; exit 1; }
mkdir -p search/dbcan logs
echo "=== Gate 1 (CAZyme) ==="
echo "  DB: $DB"
echo "  단백질: $(grep -c '^>' proteins/all_proteins.faa) 서열"
echo "  CPU: $(nproc)"
echo "  시작: $(date -Is)"
hmmsearch --domtblout search/dbcan/cazyme.domtbl --cpu "$(nproc)" -E 1e-15 --domE 1e-15 "$DB" proteins/all_proteins.faa > /dev/null
echo "  검색 완료: $(date -Is)"
python - <<'PY'
import collections, re
gh=collections.defaultdict(set); allfam=collections.defaultdict(set); n=0
for line in open('search/dbcan/cazyme.domtbl'):
    if line.startswith('#'): continue
    p=line.split()
    if len(p)<23: continue
    cov=(int(p[16])-int(p[15])+1)/max(int(p[5]),1)
    if cov<0.35: continue
    acc=p[0].split('|')[0]; fam=re.sub(r'\.hmm$','',p[3]).split('_')[0]
    allfam[acc].add(fam); n+=1
    if fam.startswith('GH'): gh[acc].add(fam)
with open('search/dbcan/overview.txt','w') as f:
    f.write('accession\tGH_families\tall_families\n')
    for acc in sorted(allfam):
        f.write(acc+'\t'+';'.join(sorted(gh.get(acc,[])))+'\t'+';'.join(sorted(allfam[acc]))+'\n')
print('  도메인 히트',n,'· 유전체',len(allfam))
PY
{ echo ""; echo "[gate1_cazyme] $(date -Is)"; echo "  DB $DB"; } >> logs/run_log.txt
echo "=== 완료 ==="
echo "다음:  python 06_call_gates.py && python 07_summarize.py"
