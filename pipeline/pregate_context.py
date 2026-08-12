import os, re, csv, collections
WINDOW = 5000
MOBILE = re.compile(r'transposase|integrase|recombinase|IS\d|insertion sequence|'
                    r'conjugat|relaxase|mobilization|plasmid replication|repA|'
                    r'type IV secretion|traA|traB|traE|excisionase', re.I)
print("=== 5단계: 이동성 요소 인접 판정 ===")
hits = collections.defaultdict(dict)
for db, path in [('CARD','search/pregate/CARD_hits.tsv'), ('VFDB','search/pregate/VFDB_hits.tsv')]:
    if not os.path.exists(path):
        print("  !", path, "없음"); continue
    n = 0
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 8: continue
        acc, prot = p[0].split('|', 1)
        hits[acc].setdefault(prot, []).append((db, p[1])); n += 1
    print("  %s: %d 히트" % (db, n))
print("  히트를 가진 유전체 %d 개" % len(hits))
print("  GFF 위치 색인")
gff_of = {}
for root, dirs, files in os.walk('genomes'):
    if 'genomic.gff' in files:
        m = re.search(r'GC[AF]_[0-9]+\.[0-9]+', root)
        if m: gff_of[m.group(0)] = os.path.join(root, 'genomic.gff')
print("    %d 개 GFF" % len(gff_of))
rows = []
summary = collections.defaultdict(lambda: {'hits':0, 'mobile':0})
for acc in sorted(hits):
    gf = gff_of.get(acc)
    if not gf:
        for prot in hits[acc]:
            summary[acc]['hits'] += 1
            rows.append([acc, prot, ';'.join(d for d,_ in hits[acc][prot]), '', '', 'no GFF', ''])
        continue
    feats = []
    for line in open(gf):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS': continue
        pid = ''; prod = ''
        m = re.search(r'protein_id=([^;]+)', f[8])
        if m: pid = m.group(1)
        m = re.search(r'product=([^;]+)', f[8])
        if m: prod = m.group(1)
        feats.append((f[0], int(f[3]), int(f[4]), prod, pid))
    by_contig = collections.defaultdict(list)
    for c,s,e,prod,pid in feats: by_contig[c].append((s,e,prod,pid))
    for c in by_contig: by_contig[c].sort()
    loc = dict((pid,(c,s,e)) for c,s,e,prod,pid in feats if pid)
    for prot, labels in hits[acc].items():
        summary[acc]['hits'] += 1
        if prot not in loc:
            rows.append([acc, prot, ';'.join(d for d,_ in labels), '', '', 'not located', '']); continue
        c,s,e = loc[prot]
        near = [prod for (ss,ee,prod,_) in by_contig[c]
                if ee >= s-WINDOW and ss <= e+WINDOW and MOBILE.search(prod or '')]
        if near:
            summary[acc]['mobile'] += 1
            rows.append([acc, prot, ';'.join(d for d,_ in labels), c, "%d-%d"%(s,e),
                         'MOBILE-ADJACENT', '; '.join(sorted(set(near))[:3])])
        else:
            rows.append([acc, prot, ';'.join(d for d,_ in labels), c, "%d-%d"%(s,e),
                         'no mobile element within 10 kb', ''])
os.makedirs('results', exist_ok=True)
w = csv.writer(open('results/pregate_context.tsv','w',newline=''), delimiter='\t')
w.writerow(['assembly','protein','database','contig','coordinates','context_call','adjacent_mobile_features'])
w.writerows(rows)
w2 = csv.writer(open('results/pregate_per_genome.tsv','w',newline=''), delimiter='\t')
w2.writerow(['assembly','determinants_detected','mobile_adjacent','pre_gate_call'])
fail = 0
for acc in sorted(summary):
    d = summary[acc]
    call = 'FAIL (transferable)' if d['mobile'] else 'pass (intrinsic only)'
    if d['mobile']: fail += 1
    w2.writerow([acc, d['hits'], d['mobile'], call])
print()
print("  유전체 %d · 전달 가능 판정 %d · 통과 %d" % (len(summary), fail, len(summary)-fail))
print("  → results/pregate_context.tsv")
print("  → results/pregate_per_genome.tsv")
open('logs/run_log.txt','a').write("\n[pregate_context] window 10000 bp · %d/%d genomes mobile-adjacent\n" % (fail, len(summary)))
