#!/usr/bin/env python3
# ============================================================
# 01c_fix_decoys.py  —  애매한 참조서열 걸러내고 다양성 확인
# ============================================================
# 실행:  python 01c_fix_decoys.py
#
# 왜 필요한가
#   NCBI 는 group II 탈탄산효소를 흔히 이렇게만 적어 둡니다.
#     "aminotransferase class I/II-fold pyridoxal phosphate-dependent enzyme"
#   이 이름만으로는 라이신인지 타이로신인지 알 수 없습니다.
#   decoy 자리에 이런 서열이 들어가면 진짜 TyrDC 가 decoy 로 잡혀
#   양성 균주를 음성으로 판정하게 됩니다.
#
#   그래서 이름이 분명하지 않은 서열은 decoy 에서 빼고,
#   대신 기능이 확립된 대장균·유산균 참조로 채웁니다.
# ============================================================

import os
import re
import time
import datetime
import requests
from collections import Counter, defaultdict

FAA = "refs/gate2_reference.faa"
BAK = "refs/gate2_reference.faa.bak"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# decoy 에 들어가면 안 되는 애매한 이름
VAGUE = re.compile(
    r"aminotransferase class|pyridoxal phosphate-dependent enzyme$|"
    r"hypothetical|uncharacteri|unnamed|DUF\d+", re.I)

# 계열마다 이름에 반드시 있어야 할 낱말
REQUIRED = {
    "TARGET_tyrDC": ["tyrosine decarboxylase", "tyrdc", "tdca"],
    "TARGET_tyrP":  ["tyrosine permease", "tyrp"],
    "TARGET_aadc":  ["aromatic", "decarboxylase"],
    "DECOY_gadB":   ["glutamate decarboxylase", "gadb"],
    "DECOY_ldcA":   ["lysine decarboxylase", "ldca", "cada", "ldcc"],
    "DECOY_odcA":   ["ornithine decarboxylase", "odca", "spef"],
    "PREGATE_hdcA": ["histidine decarboxylase", "hdca"],
}

# 기능이 확립된 대체 참조 — 이름이 분명한 것만 씁니다
REPLACE = {
    "DECOY_ldcA": [
        'lysine decarboxylase, constitutive AND Escherichia coli[Organism]',
        'lysine decarboxylase, inducible AND Escherichia coli[Organism]',
        '"lysine decarboxylase" AND Lactobacillus[Organism] NOT hypothetical',
        '"lysine decarboxylase" AND Enterococcus[Organism] NOT hypothetical',
    ],
    "DECOY_odcA": [
        'ornithine decarboxylase, inducible AND Escherichia coli[Organism]',
        '"ornithine decarboxylase" AND Lactobacillus[Organism] NOT hypothetical',
    ],
    "DECOY_gadB": [
        'glutamate decarboxylase AND Escherichia coli[Organism]',
        '"glutamate decarboxylase" AND Lactobacillus[Organism] NOT hypothetical',
    ],
    "TARGET_tyrP": [
        '"tyrosine permease" AND Enterococcus[Organism] NOT hypothetical',
        '"tyrosine permease" AND Lactobacillus[Organism] NOT hypothetical',
    ],
}


def read_fasta(path):
    recs, head, seq = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if head:
                recs.append((head, "".join(seq)))
            head, seq = line.strip()[1:], []
        else:
            seq.append(line.strip())
    if head:
        recs.append((head, "".join(seq)))
    return recs


def esearch(term, n):
    r = requests.get(f"{EUTILS}/esearch.fcgi", timeout=60, params={
        "db": "protein", "term": term, "retmax": n, "retmode": "json"})
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def efetch(ids):
    if not ids:
        return ""
    r = requests.get(f"{EUTILS}/efetch.fcgi", timeout=120, params={
        "db": "protein", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"})
    r.raise_for_status()
    return r.text


def main():
    if not os.path.exists(FAA):
        raise SystemExit(f"파일이 없습니다: {FAA}")

    recs = read_fasta(FAA)
    if not os.path.exists(BAK):
        with open(BAK, "w") as f:
            for h, s in recs:
                f.write(f">{h}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i+60] + "\n")
        print(f"원본을 백업했습니다: {BAK}\n")

    # ------------------------------------------------------------
    # 1) 이름이 애매하거나 계열과 맞지 않는 서열을 걸러냅니다
    # ------------------------------------------------------------
    keep, dropped = [], defaultdict(list)
    for head, seq in recs:
        label = head.split("|")[0]
        desc = " ".join(head.split()[1:])
        req = REQUIRED.get(label, [])
        ok_name = any(k in desc.lower() for k in req) if req else True
        vague = bool(VAGUE.search(desc))
        if vague or not ok_name:
            dropped[label].append((head.split("|")[1].split()[0], desc[:56]))
        else:
            keep.append((head, seq))

    print("=== 1c 단계: 참조서열 정리 ===")
    print("걸러낸 서열")
    if not dropped:
        print("  없음 — 모두 이름이 분명합니다")
    for lab, items in sorted(dropped.items()):
        print(f"  {lab}  {len(items)}건")
        for acc, d in items[:4]:
            print(f"      {acc:18s} {d}")

    # ------------------------------------------------------------
    # 2) 부족해진 계열을 이름이 분명한 참조로 채웁니다
    # ------------------------------------------------------------
    have = Counter(h.split("|")[0] for h, _ in keep)
    seen = {h.split("|")[1].split()[0] for h, _ in keep}
    added = Counter()

    for lab, terms in REPLACE.items():
        need = 5 - have[lab]
        if need <= 0:
            continue
        for term in terms:
            if added[lab] >= need:
                break
            try:
                ids = esearch(term, 10)
                time.sleep(0.34)
                fasta = efetch(ids)
                time.sleep(0.34)
            except Exception as e:
                print(f"  ! {lab}: {e}")
                continue
            for block in fasta.split("\n>"):
                if added[lab] >= need:
                    break
                block = block.strip()
                if not block:
                    continue
                if not block.startswith(">"):
                    block = ">" + block
                lines = block.split("\n")
                h0, s0 = lines[0], "".join(lines[1:])
                acc = h0[1:].split()[0]
                desc = " ".join(h0[1:].split()[1:])
                if acc in seen or len(s0) < 200:
                    continue
                if VAGUE.search(desc):
                    continue
                if not any(k in desc.lower() for k in REQUIRED.get(lab, [])):
                    continue
                seen.add(acc)
                keep.append((f"{lab}|{acc} {desc}", s0))
                added[lab] += 1

    if added:
        print("\n보충한 서열")
        for lab, n in sorted(added.items()):
            print(f"  {lab:16s} +{n}")

    # ------------------------------------------------------------
    # 3) 저장하고 다양성을 확인합니다
    # ------------------------------------------------------------
    with open(FAA, "w") as f:
        for h, s in keep:
            f.write(f">{h}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")

    final = Counter(h.split("|")[0] for h, _ in keep)
    print("\n최종 구성")
    warn = []
    for lab in sorted(REQUIRED):
        n = final[lab]
        mark = "  " if n >= 3 else "★ "
        print(f"  {mark}{lab:16s} {n:3d} 건")
        if n < 2:
            warn.append(lab)

    # 표적 계열의 균종 다양성 — 한 연구에서 나온 중복 서열이면 다양성이 없습니다
    print("\nTARGET_tyrDC 의 균종 분포")
    orgs = Counter()
    for h, _ in keep:
        if h.startswith("TARGET_tyrDC"):
            m = re.search(r"\[([^\]]+)\]", h)
            if m:
                orgs[m.group(1)] += 1
    for o, n in orgs.most_common():
        print(f"    {o:44s} {n}")
    if len(orgs) < 2:
        print("    ★ 균종이 하나뿐입니다. 다양성이 부족하면 다른 속의 TyrDC 를 놓칠 수 있습니다.")

    with open("logs/run_log.txt", "a") as lg:
        lg.write(f"\n[01c_fix_decoys] {datetime.datetime.now().isoformat()}\n")
        for lab, items in dropped.items():
            lg.write(f"  제거 {lab}: {len(items)}\n")
        for lab, n in added.items():
            lg.write(f"  보충 {lab}: +{n}\n")
        lg.write(f"  최종 {dict(final)}\n")

    print()
    if warn:
        print(f"★ 아직 부족한 계열: {', '.join(warn)}")
        print("  이 상태로 진행하려면 06_call_gates.py 의 MARGIN 을 10 에서 20 으로 올리십시오.")
    else:
        print("모든 계열이 충분합니다.")
    print("\n다음:  python 02_build_panel.py --taxa taxa_list.txt --out panel/")


if __name__ == "__main__":
    main()
