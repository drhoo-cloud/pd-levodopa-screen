#!/usr/bin/env bash
# ============================================================
# gate1_cazyme.sh  —  Gate 1 (CAZyme) 주석과 GH family 집계
# ============================================================
# 실행:  bash gate1_cazyme.sh
#
# ------------------------------------------------------------
# 구 버전에서 바뀐 점  ★ 이것이 이 파일을 다시 쓴 이유입니다
# ------------------------------------------------------------
#   구 버전은 dbCAN.hmm 에 hmmsearch 를 직접 걸었습니다.
#   즉 dbCAN 의 세 가지 예측 방법 중 HMMER 하나만 쓴 것입니다.
#
#   dbCAN 은 단백질마다 세 가지로 독립 판정합니다.
#       HMMER      dbCAN-HMM 프로파일 검색
#       DIAMOND    CAZy 서열 데이터베이스 검색
#       dbCAN_sub  하위 계열 프로파일
#
#   한 가지만 쓰면 도메인 히트가 과다 계상됩니다.
#   실제로 L. casei 33개 유전체의 GH family 중앙값이 26 으로 나왔는데,
#   세 방법 중 둘 이상이 일치하는 것만 세면 21 입니다.
#   종별 중앙값이 몇 family 씩 밀립니다.
#
#   그래서 이 스크립트는 run_dbcan 을 --tools all 로 돌리고,
#   overview.txt 의 '#ofTools' 가 2 이상인 행만 집계합니다.
# ------------------------------------------------------------
set -euo pipefail

PROT=proteins/all_proteins.faa
OUTDIR=search/dbcan
MIN_TOOLS=2                     # 세 방법 중 몇 개가 일치해야 인정할지

[ -f "$PROT" ] || { echo "단백질 파일 없음: $PROT"; exit 1; }
[ -n "${DBCAN_DB:-}" ] || { echo "DBCAN_DB 가 설정되지 않았습니다."; echo "  export DBCAN_DB=\$PWD/db"; exit 1; }

if ! command -v run_dbcan >/dev/null 2>&1; then
  cat <<'MSG'
run_dbcan 이 없습니다. 아래로 설치하십시오.

    pip install dbcan
    mkdir -p db && cd db
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt
    hmmpress dbCAN-HMMdb-V12.txt
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa
    diamond makedb --in CAZyDB.07262023.fa -d CAZy
    wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN_sub.hmm
    hmmpress dbCAN_sub.hmm
    cd ..
    export DBCAN_DB=$PWD/db

DIAMOND 와 dbCAN_sub 데이터베이스가 모두 있어야 세 방법이 돌아갑니다.
HMMER 것만 받으면 구 버전과 같은 결과가 나옵니다.
MSG
  exit 1
fi

mkdir -p "$OUTDIR" logs
echo "=== Gate 1 (CAZyme) ==="
echo "  단백질: $(grep -c '^>' "$PROT") 서열"
echo "  DB:     $DBCAN_DB"
echo "  방법:   HMMER + DIAMOND + dbCAN_sub (일치 $MIN_TOOLS 개 이상만 인정)"
echo "  CPU:    $(nproc)"
echo "  시작:   $(date -Is)"

run_dbcan "$PROT" protein \
  --out_dir "$OUTDIR" \
  --db_dir "$DBCAN_DB" \
  --tools all \
  --hmm_cpu "$(nproc)" \
  --dia_cpu "$(nproc)" \
  --dbcan_thread "$(nproc)"

echo "  검색 완료: $(date -Is)"

OV="$OUTDIR/overview.txt"
[ -f "$OV" ] || { echo "overview.txt 가 없습니다: $OV"; exit 1; }

# ------------------------------------------------------------
# 집계 — overview.txt 를 열 이름으로 읽습니다.
#        dbCAN 버전에 따라 열 순서가 바뀌므로 위치로 읽으면 안 됩니다.
# ------------------------------------------------------------
MIN_TOOLS="$MIN_TOOLS" OUTDIR="$OUTDIR" python3 - <<'PY'
import os, re, csv, collections

OUTDIR = os.environ["OUTDIR"]
MIN_TOOLS = int(os.environ["MIN_TOOLS"])
OV = os.path.join(OUTDIR, "overview.txt")

def norm(s):
    return re.sub(r"[^a-z]", "", s.lower())

FAM = re.compile(r"^[A-Z]{2,3}\d+$")

def parse_fams(cell):
    """\uc140 \ud558\ub098\uc5d0\uc11c family \uc774\ub984\uc744 \ubf51\uc2b5\ub2c8\ub2e4.

    'GH73_e145|CBM50_e1043'  -> {GH73, CBM50}
    'GH13_20(176-480)'       -> {GH13}
    \uad6c\ubd84\uc790\uc5d0 \ud30c\uc774\ud504(|)\uac00 \ubc18\ub4dc\uc2dc \ub4e4\uc5b4\uac00\uc57c \ud569\ub2c8\ub2e4.
    Recommend Results \uc5f4\uc774 \uadf8\uac83\uc744 \uad6c\ubd84\uc790\ub85c \uc501\ub2c8\ub2e4.
    """
    out = set()
    if not cell or str(cell).strip() in ("-", "N", "NA"):
        return out
    for tok in re.split(r"[+;,|\s]+", str(cell)):
        tok = tok.split("(")[0].split("_")[0].strip()
        if FAM.fullmatch(tok):
            out.add(tok)
    return out

with open(OV) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    idx = {norm(h): i for i, h in enumerate(header)}

    # 열 이름은 dbCAN 버전마다 조금씩 다릅니다
    def find(*cands):
        for c in cands:
            if norm(c) in idx:
                return idx[norm(c)]
        return None

    i_gene  = find("Gene ID", "GeneID")
    i_tools = find("#ofTools", "Number of Tools", "#ofTool")
    i_rec   = find("Recommend Results", "RecommendResults")
    i_hmm   = find("dbCAN_hmm", "HMMER", "dbCAN")
    i_sub   = find("dbCAN_sub", "eCAMI", "dbCANsub")
    i_dia   = find("DIAMOND", "Diamond")

    if i_gene is None:
        raise SystemExit(f"Gene ID \uc5f4\uc744 \ucc3e\uc9c0 \ubabb\ud588\uc2b5\ub2c8\ub2e4. \ud5e4\ub354: {header}")
    if i_tools is None:
        print("  ! #ofTools \uc5f4\uc774 \uc5c6\uc2b5\ub2c8\ub2e4. \uc138 \ubc29\ubc95 \uc5f4\uc5d0\uc11c \uc9c1\uc811 \uc149\ub2c8\ub2e4.")
    if i_rec is None:
        print("  ! Recommend Results \uc5f4\uc774 \uc5c6\uc2b5\ub2c8\ub2e4. \uc138 \ubc29\ubc95\uc758 \ud569\uc9d1\ud569\uc744 \uc501\ub2c8\ub2e4.")
        print("    \uc774 \uacbd\uc6b0 family \uac1c\uc218\uac00 \ub192\uac8c \ub098\uc635\ub2c8\ub2e4. dbCAN \ubc84\uc804\uc744 \ud655\uc778\ud558\uc2ed\uc2dc\uc624.")

    gh      = collections.defaultdict(set)
    allfam  = collections.defaultdict(set)
    dropped = 0
    kept    = 0
    fallback = 0

    for row in reader:
        if not row or len(row) <= i_gene:
            continue

        # --- \uba87 \uac1c \ubc29\ubc95\uc774 \uc774 \ub2e8\ubc31\uc9c8\uc744 CAZyme \uc73c\ub85c \ubd88\ub800\ub294\uac00 -------------
        if i_tools is not None:
            try:
                ntools = int(str(row[i_tools]).strip())
            except ValueError:
                continue
        else:
            ntools = 0
            for i in (i_hmm, i_sub, i_dia):
                if i is not None and i < len(row):
                    v = str(row[i]).strip()
                    if v and v not in ("-", "N", "NA"):
                        ntools += 1

        if ntools < MIN_TOOLS:
            dropped += 1
            continue
        kept += 1

        # --- family \uc774\ub984 \ubf51\uae30 -----------------------------------------
        #   \uc138 \ubc29\ubc95\uc774 \uac19\uc740 \ub2e8\ubc31\uc9c8\uc744 \uc11c\ub85c \ub2e4\ub978 family \ub85c \ubd80\ub97c \uc218 \uc788\uc2b5\ub2c8\ub2e4.
        #   \uc608:  dbCAN_hmm=GH179  dbCAN_sub=GH109  DIAMOND=-   (#ofTools=2)
        #   \ud569\uc9d1\ud569\uc744 \uc138\uba74 \ud55c \ub2e8\ubc31\uc9c8\uc774 \ub450 family \ub85c \uacc4\uc0b0\ub418\uc5b4 \uac1c\uc218\uac00 \ubd80\ud480\ub824\uc9d1\ub2c8\ub2e4.
        #   dbCAN \uc740 'Recommend Results' \uc5f4\uc5d0 \uc790\uccb4 \uc870\uc815 \uacb0\uacfc\ub97c \ub0b4\ub193\uc2b5\ub2c8\ub2e4.
        #   \uadf8\uac83\uc744 \ub530\ub985\ub2c8\ub2e4. \uad6c\ubd84\uc790\ub294 \ud30c\uc774\ud504(|) \uc785\ub2c8\ub2e4.
        fams = set()
        if i_rec is not None and i_rec < len(row):
            fams = parse_fams(row[i_rec])
        if not fams:
            fallback += 1
            for i in (i_hmm, i_sub, i_dia):
                if i is not None and i < len(row):
                    fams |= parse_fams(row[i])

        fams.discard("GH0")     # dbCAN \uc758 \ubbf8\ubd84\ub958 \uc790\ub9ac\ud45c\uc2dc\uc790. \uc2e4\uc81c family \uac00 \uc544\ub2d9\ub2c8\ub2e4.
        if not fams:
            continue
        acc = row[i_gene].split("|")[0]
        allfam[acc] |= fams
        gh[acc] |= {f for f in fams if f.startswith("GH")}

out = os.path.join(OUTDIR, "overview_consensus.txt")
with open(out, "w") as f:
    f.write("accession\tn_GH_families\tGH_families\tall_families\n")
    for acc in sorted(allfam):
        g = sorted(gh.get(acc, []))
        f.write(f"{acc}\t{len(g)}\t{';'.join(g)}\t{';'.join(sorted(allfam[acc]))}\n")

print(f"  일치 {MIN_TOOLS}개 이상: {kept} 단백질 · 버림 {dropped} 단백질")
print(f"  유전체 {len(allfam)} 개 → {out}")

if gh:
    import statistics
    counts = [len(v) for v in gh.values()]
    print(f"  GH family 수: 중앙값 {statistics.median(counts):g} "
          f"· 최소 {min(counts)} · 최대 {max(counts)}")
PY

{
  echo ""
  echo "[gate1_cazyme] $(date -Is)"
  echo "  run_dbcan --tools all · DB $DBCAN_DB"
  echo "  일치 기준 #ofTools >= $MIN_TOOLS"
  echo "  출력 $OUTDIR/overview_consensus.txt"
} >> logs/run_log.txt

echo "=== Gate 1 완료 ==="
echo "  구 버전 결과와 반드시 비교하십시오. 중앙값이 몇 family 낮아야 정상입니다."
echo "다음:  python 06_call_gates.py && python 07_summarize.py"
