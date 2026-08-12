import os, re, time
t0 = time.time()
print("1) protein.faa 파일 한 번에 찾는 중")
idx = {}
for root, dirs, files in os.walk("genomes"):
    if "protein.faa" in files:
        m = re.search(r'GC[AF]_[0-9]+\.[0-9]+', root)
        if m:
            idx[m.group(0)] = os.path.join(root, "protein.faa")
print(f"   {len(idx)}개 발견  ({time.time()-t0:.0f}초)")

want = [l.strip() for l in open("panel/_want.txt") if l.strip()]
print(f"2) 병합 중 (요청 {len(want)}건)")
found = missing = nseq = 0
with open("proteins/all_proteins.faa", "w") as out, \
     open("logs/missing_protein_faa.txt", "w") as miss:
    for i, acc in enumerate(want, 1):
        p = idx.get(acc)
        if not p:
            miss.write(acc + "\n"); missing += 1; continue
        with open(p) as f:
            for line in f:
                if line.startswith(">"):
                    out.write(">" + acc + "|" + line[1:]); nseq += 1
                else:
                    out.write(line)
        found += 1
        if i % 200 == 0:
            print(f"   {i}/{len(want)} ...")

print()
print("=== 병합 완료 ===")
print(f"  확보 {found} / {len(want)} · 누락 {missing}")
print(f"  단백질 {nseq} 서열 → proteins/all_proteins.faa")
print(f"  걸린 시간 {time.time()-t0:.0f}초")
with open("logs/run_log.txt", "a") as lg:
    lg.write(f"\n[fast_merge] 확보 {found}/{len(want)} 누락 {missing} 서열 {nseq}\n")
print("\n다음:  python 04_validate_controls.py")
