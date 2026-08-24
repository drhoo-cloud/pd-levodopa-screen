#!/usr/bin/env bash
# ============================================================
# pregate.sh  —  프리게이트: 안전성 선별
# ============================================================
# 실행:  bash pregate.sh [accession_목록파일]
#        기본값 panel/accessions_ent.txt
#
# 무엇을 하나
#   (A) 전달 가능한 내성·병독 인자   CARD / VFDB 검색
#   (B) 히스티딘 탈탄산효소           바이오제닉 아민 안전성
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점  ★ (B) 가 통째로 빠져 있었습니다
# ------------------------------------------------------------
#   원고 2절은 프리게이트가 히스티딘 탈탄산효소도 본다고 기술합니다.
#   그런데 구 pregate.sh 는 CARD 와 VFDB 만 검색했습니다.
#   HDC 참조 서열은 01c_fix_decoys.py 안에만 있었고 어디서도 쓰이지 않았습니다.
#   논문과 코드가 어긋나 있던 지점이라 여기서 메웁니다.
#
#   왜 Gate 2 가 아니라 프리게이트인가
#     젖산균의 HDC 는 pyruvoyl 의존형입니다. group II 피리독살인산 계열이 아닙니다.
#     Gate 2 의 프로파일로는 잡히지 않습니다. 따라서 여기서 따로 잡아야 합니다.
#
#   왜 Gate 2 와 다른 방식(상동성 검색)을 쓰는가
#     Gate 2 는 "levodopa 를 가로채는가" 라는 이분 판정입니다.
#     HDC 선별은 "히스타민 생성 가능성이 있는가" 라는 안전 경고입니다.
#     경고는 놓치는 쪽이 잘못 울리는 쪽보다 나쁘므로 임계를 느슨하게 둡니다.
#     이 결과는 판정이 아니라 회부 사유입니다.
# ------------------------------------------------------------
set -uo pipefail

LIST="${1:-panel/accessions_ent.txt}"
CARD=db/amr/protein_fasta_protein_homolog_model.fasta
VFDB=db/amr/vfdb.fas
HDC=refs/pregate_hdc.faa

# 임계값 — 사전 확정
AMR_ID=80;  AMR_COV=70;  AMR_E=1e-10    # 내성·병독: 엄격 (같은 유전자여야 함)
HDC_ID=40;  HDC_COV=70;  HDC_E=1e-10    # HDC: 느슨 (계열 수준 경고)

[ -f "$LIST" ] || { echo "accession 목록이 없습니다: $LIST"; exit 1; }
mkdir -p search/pregate results logs

N=$(grep -c . "$LIST")
echo "=== 프리게이트: 안전성 선별 ==="
echo "  대상 $N 균주   ($LIST)"

# ------------------------------------------------------------
# 대상 유전체의 단백질만 뽑아냅니다
# ------------------------------------------------------------
echo "  단백질 추출"
LIST="$LIST" python3 - <<'PY'
import os
keep = set(l.strip() for l in open(os.environ['LIST']) if l.strip())
out = open('search/pregate/ent.faa', 'w')
write, n = False, 0
for line in open('proteins/all_proteins.faa'):
    if line[0] == '>':
        write = line[1:].split('|')[0] in keep
        if write:
            n += 1
    if write:
        out.write(line)
out.close()
print(f'    {n} 서열')
PY

# ============================================================
# (A) 전달 가능한 내성·병독 인자
# ============================================================
echo
echo "  [A] 내성·병독 인자"
for DB in CARD VFDB; do
  [ "$DB" = CARD ] && F=$CARD || F=$VFDB
  if [ ! -f "$F" ]; then
    echo "    ! $DB 데이터베이스 없음: $F — 건너뜁니다"
    continue
  fi
  echo "    $DB 검색"
  diamond makedb --in "$F" -d "search/pregate/${DB}" --quiet
  diamond blastp -q search/pregate/ent.faa -d "search/pregate/${DB}" \
    -o "search/pregate/${DB}_hits.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --evalue "$AMR_E" --id "$AMR_ID" --query-cover "$AMR_COV" \
    --max-target-seqs 5 --threads "$(nproc)" --quiet
  echo "      히트 $(wc -l < "search/pregate/${DB}_hits.tsv") 건"
done

echo
echo "    ※ 히트가 있다는 것만으로는 탈락이 아닙니다."
echo "       이동성 요소에 인접해 있어야 '전달 가능' 입니다."
echo "       판정은 pregate_context.py 가 합니다."

# ============================================================
# (B) 히스티딘 탈탄산효소  ★ 새로 추가된 블록
# ============================================================
echo
echo "  [B] 히스티딘 탈탄산효소 (바이오제닉 아민)"
if [ ! -f "$HDC" ]; then
  echo "    ! 참조 서열 없음: $HDC"
  echo "      python 01_reference_proteins.py 를 먼저 실행하십시오."
  echo "      (UniProt 에서 reviewed 세균 HDC 를 받아 이 파일을 만듭니다)"
else
  NREF=$(grep -c '^>' "$HDC")
  echo "    참조 $NREF 서열 · 임계 identity>=${HDC_ID}% coverage>=${HDC_COV}%"
  diamond makedb --in "$HDC" -d search/pregate/HDC --quiet
  diamond blastp -q search/pregate/ent.faa -d search/pregate/HDC \
    -o search/pregate/HDC_hits.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --evalue "$HDC_E" --id "$HDC_ID" --query-cover "$HDC_COV" \
    --max-target-seqs 5 --threads "$(nproc)" --quiet
  echo "      히트 $(wc -l < search/pregate/HDC_hits.tsv) 건"

  # 유전체 단위로 정리 — 가장 강한 히트 하나만 남깁니다
  LIST="$LIST" python3 - <<'PY'
import csv, collections, os

path = 'search/pregate/HDC_hits.tsv'
best = {}
for line in open(path):
    p = line.rstrip('\n').split('\t')
    if len(p) < 8:
        continue
    acc = p[0].split('|')[0]
    rec = {'protein': p[0], 'subject': p[1], 'pident': float(p[2]),
           'evalue': p[6], 'bitscore': float(p[7])}
    if acc not in best or rec['bitscore'] > best[acc]['bitscore']:
        best[acc] = rec

# 대상 유전체 전체를 분모에 넣습니다 (히트 없음 = 명백한 음성)
keep = [l.strip() for l in open(os.environ.get('LIST', 'panel/accessions_ent.txt')) if l.strip()]

os.makedirs('results', exist_ok=True)
with open('results/pregate_hdc.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['assembly', 'hdc_call', 'best_identity_pct',
                'best_subject', 'protein', 'evalue'])
    flagged = 0
    for acc in sorted(keep):
        r = best.get(acc)
        if r:
            flagged += 1
            w.writerow([acc, 'FLAG (HDC homolog)', f"{r['pident']:.1f}",
                        r['subject'], r['protein'], r['evalue']])
        else:
            w.writerow([acc, 'none detected', '', '', '', ''])

print(f"      HDC 상동체 검출 {flagged} / {len(keep)} 유전체")
print( "      → results/pregate_hdc.tsv")
if flagged:
    print()
    print("      ★ FLAG 는 탈락이 아니라 회부입니다.")
    print("        히스타민 생성은 배양 실험으로 확인해야 합니다.")
    print("        젖산균 HDC 는 pyruvoyl 의존형이라 Gate 2 로는 잡히지 않습니다.")
PY
fi

# ============================================================
{
  echo ""
  echo "[pregate] $(date -Is)"
  echo "  대상 $N 균주 ($LIST)"
  echo "  A) CARD $(wc -l < search/pregate/CARD_hits.tsv 2>/dev/null || echo 0) 건 · VFDB $(wc -l < search/pregate/VFDB_hits.tsv 2>/dev/null || echo 0) 건"
  echo "     기준 identity>=$AMR_ID coverage>=$AMR_COV evalue<=$AMR_E"
  echo "  B) HDC $(wc -l < search/pregate/HDC_hits.tsv 2>/dev/null || echo 0) 건"
  echo "     기준 identity>=$HDC_ID coverage>=$HDC_COV evalue<=$HDC_E (경고용, 느슨함)"
} | tee -a logs/run_log.txt

echo
echo "=== 프리게이트 검색 완료 ==="
echo "다음:  python pregate_context.py     (이동성 요소 인접 판정)"
