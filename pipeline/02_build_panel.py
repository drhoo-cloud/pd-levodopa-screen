#!/usr/bin/env python3
# ============================================================
# 02_build_panel.py  —  대상 유전체 목록 확정
# ============================================================
# 실행:  python 02_build_panel.py --taxa taxa_list.txt --out panel/
#
# 무엇을 하나
#   1) 종 단위로 NCBI RefSeq 어셈블리 전체를 받습니다.
#   2) 균주명을 정규화해 우리 패널 목록과 조인합니다.
#   3) 결과를 네 상태 중 하나로 분류합니다.
#
#        RESOLVED            accession 확보 — 분석 대상
#        AMBIGUOUS_MULTIPLE  후보가 둘 이상 — 사람이 골라야 함
#        NO_PUBLIC_GENOME    균주명은 있으나 공개 유전체 없음
#        NOT_NAMED           원문에 균주 표기 없음 — 해소 불가, 이것이 결과
#
# ★ 이 스크립트는 절대 하지 않는 것
#   · 후보가 여럿일 때 임의로 하나를 고르지 않습니다.
#   · 조회 실패 시 같은 종의 다른 어셈블리로 대체하지 않습니다.
#   위 둘을 하는 순간 "균주 단위 판정" 이라는 논문의 전제가 무너집니다.
# ============================================================

import os
import re
import csv
import json
import time
import argparse
import datetime
import requests

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
TODAY = datetime.date.today().isoformat()

# 같은 균주의 다른 표기 — 조인 실패 시 이 값으로 재시도합니다
ALIASES = {
    "gg": ["atcc53103", "lgg"],
    "atcc53103": ["gg", "lgg"],
    "299v": ["dsm9843"],
    "dsm9843": ["299v"],
    "lp01": ["lmgp21021"],
    "ls01": ["dsm22775"],
    "la02": ["dsm21717"],
    "lr06": ["dsm21981"],
    "atcc14917": ["dsm20174", "jcm1149", "lmg6907", "cgmcc12437", "ncimb11974", "lp39"],
    "ps128": ["dsm28632"],
    "shirota": ["yit9029", "lcs"],
    "wcfs1": ["ncimb8826"],
}


def norm(s):
    """비교용 정규화: 소문자 + 괄호 안 제거 + 공백/하이픈/점 제거"""
    if not s:
        return ""
    s = str(s).lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(r"[ \-_.=]", "", s)
    return s.strip()


def fetch_taxon(taxon, page_size=1000):
    """한 종의 RefSeq 어셈블리를 전부 받습니다 (페이지 넘김 포함)"""
    out, token = [], None
    while True:
        params = {"filters.assembly_source": "refseq",
                  "filters.exclude_atypical": "true",
                  "page_size": page_size}
        if token:
            params["page_token"] = token
        url = f"{API}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report"
        r = requests.get(url, params=params, timeout=90)
        if r.status_code != 200:
            print(f"    ! HTTP {r.status_code}")
            break
        js = r.json()
        for rep in js.get("reports", []):
            org = rep.get("organism", {}) or {}
            info = rep.get("assembly_info", {}) or {}
            out.append({
                "accession": rep.get("accession", ""),
                "organism": org.get("organism_name", ""),
                "strain": (org.get("infraspecific_names", {}) or {}).get("strain", ""),
                "level": info.get("assembly_level", ""),
                "biosample": (info.get("biosample", {}) or {}).get("accession", ""),
                "isolation_source": "",
            })
        token = js.get("next_page_token")
        if not token:
            break
        time.sleep(0.34)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxa", required=True, help="종 목록 파일 (한 줄에 하나)")
    ap.add_argument("--strains", default="", help="특정 균주 목록 CSV (선택). 열: strain")
    ap.add_argument("--out", default="panel/")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    taxa = [l.strip() for l in open(args.taxa) if l.strip() and not l.startswith("#")]
    print(f"=== 2단계: 유전체 목록 확보 ({len(taxa)} 종) ===\n")

    all_recs, index = [], {}
    for t in taxa:
        recs = fetch_taxon(t)
        print(f"  {t:38s} {len(recs):6d} assemblies")
        all_recs.extend(recs)
        for rec in recs:
            for k in {norm(rec["strain"]), norm(rec["organism"].split()[-1])}:
                if k:
                    index.setdefault(k, []).append(rec)
        time.sleep(0.34)

    # 전체 목록 저장 — 이것이 종 단위 패널의 원자료가 됩니다
    all_tsv = os.path.join(args.out, "assemblies_all.tsv")
    with open(all_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_recs[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(all_recs)
    print(f"\n  전체 {len(all_recs)}건 → {all_tsv}")

    # ------------------------------------------------------------
    # 특정 균주 목록이 주어지면 조인해서 상태를 매깁니다
    # ------------------------------------------------------------
    if args.strains and os.path.exists(args.strains):
        rows, counts = [], {}
        with open(args.strains, newline="") as f:
            for row in csv.DictReader(f):
                st = (row.get("strain") or "").strip()
                if not st:
                    continue
                if "not stated" in st.lower():
                    status, acc, note = "NOT_NAMED", "", "원문 미표기 — 감사의 분모"
                else:
                    key = norm(st)
                    cands = list(index.get(key, []))
                    if not cands:
                        for alt in ALIASES.get(key, []):
                            cands.extend(index.get(alt, []))
                    seen, uniq = set(), []
                    for c in cands:
                        if c["accession"] not in seen:
                            seen.add(c["accession"]); uniq.append(c)
                    if len(uniq) == 1:
                        status, acc = "RESOLVED", uniq[0]["accession"]
                        note = f"{uniq[0]['level']} | strain={uniq[0]['strain']} | {TODAY}"
                    elif len(uniq) > 1:
                        status, acc = "AMBIGUOUS_MULTIPLE", ""
                        note = "후보 " + str(len(uniq)) + "건: " + \
                               "; ".join(f"{c['accession']}({c['level']})" for c in uniq[:8]) + \
                               " — 하나를 고르고 이유를 적을 것"
                    else:
                        status, acc, note = "NO_PUBLIC_GENOME", "", \
                            "RefSeq 에 없음 — 원 논문 Data availability 또는 제조사 확인"
                counts[status] = counts.get(status, 0) + 1
                rows.append({"strain": st, "assembly_accession": acc,
                             "status": status, "note": note})

        res_tsv = os.path.join(args.out, "strain_resolution.tsv")
        with open(res_tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["strain", "assembly_accession", "status", "note"],
                               delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"\n  균주 해소 결과 → {res_tsv}")
        for k, v in sorted(counts.items()):
            print(f"    {k:20s} {v}")

        # 분석 대상 accession 만 뽑기
        acc_txt = os.path.join(args.out, "accessions_all.txt")
        with open(acc_txt, "w") as f:
            for r in rows:
                if r["status"] == "RESOLVED":
                    f.write(r["assembly_accession"] + "\n")
        print(f"  분석 대상 → {acc_txt}")

    with open("logs/run_log.txt", "a") as lg:
        lg.write(f"\n[02_build_panel] {datetime.datetime.now().isoformat()}\n")
        lg.write(f"  종 {len(taxa)} · 어셈블리 {len(all_recs)}\n")

    print("\n=== 2단계 완료 ===")
    print("★ strain_resolution.tsv 에서 AMBIGUOUS_MULTIPLE 행을 사람이 정리하십시오.")
    print("   정리가 끝나면:  bash 03_download_genomes.sh panel/accessions_all.txt")


if __name__ == "__main__":
    main()
