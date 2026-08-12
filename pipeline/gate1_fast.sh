#!/usr/bin/env bash
set -uo pipefail
SRC=$PWD
W=$HOME/g1work
mkdir -p "$W"
echo "=== Gate 1 (빠른 방식) ==="

echo "1) 리눅스 내부 디스크로 복사"
[ -f "$W/prot.faa" ] || cp proteins/all_proteins.faa "$W/prot.faa"
echo "   $(du -h "$W/prot.faa" | cut -f1)"

echo "2) GH 계열 HMM 만 추출"
if [ ! -f "$W/gh.hmm" ]; then
  awk '/^HMMER3/{buf="";keep=0} {buf=buf $0 ORS}
       /^NAME/{if($2 ~ /^GH/) keep=1}
       /^\/\/$/{if(keep) printf "%s", buf; buf="";keep=0}' db/dbCAN.hmm > "$W/gh.hmm"
  hmmpress "$W/gh.hmm" >/dev/null 2>&1
fi
echo "   모델 $(grep -c '^NAME' "$W/gh.hmm") 개 (전체 875개 중 GH 만)"

echo "3) 검색 시작: $(date -Is)"
cd "$W"
hmmsearch --domtblout cazyme.domtbl --cpu "$(nproc)" -E 1e-15 --domE 1e-15 --noali \
  gh.hmm prot.faa > /dev/null
echo "   완료: $(date -Is)"

cd "$SRC"
mkdir -p search/dbcan
cp "$W/cazyme.domtbl" search/dbcan/
python - <<'PY'
import collections, re
gh=collections.defaultdict(set); n=0
for line in open('search/dbcan/cazyme.domtbl'):
    if line.startswith('#'): continue
    p=line.split()
    if len(p)<23: continue
    cov=(int(p[16])-int(p[15])+1)/max(int(p[5]),1)
    if cov<0.35: continue
    acc=p[0].split('|')[0]
    fam=re.sub(r'\.hmm$','',p[3]).split('_')[0]
    if fam.startswith('GH'): gh[acc].add(fam); n+=1
with open('search/dbcan/overview.txt','w') as f:
    f.write('accession\tGH_families\tall_families\n')
    for acc in sorted(gh):
        v=';'.join(sorted(gh[acc])); f.write(acc+'\t'+v+'\t'+v+'\n')
print('  도메인 히트',n,'· 유전체',len(gh))
PY
{ echo ""; echo "[gate1_fast] $(date -Is)"; echo "  GH-only HMM · 내부 디스크"; } >> logs/run_log.txt
echo "=== 완료 ==="
echo "다음:  python 06_call_gates.py && python 07_summarize.py"
