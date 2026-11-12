import sys, math, collections
DIST, OUT_MEMB, OUT_SUM, N_EXPECT = sys.argv[1:5]
N_EXPECT = int(N_EXPECT)
THRESH = 0.005          # Mash distance <= 0.005  ==  ANI >= 99.5%
nodes, edges = set(), []
with open(DIST) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        a = p[0].split("/")[-1].replace(".fna","")
        b = p[1].split("/")[-1].replace(".fna","")
        nodes.add(a); nodes.add(b)
        if a != b and float(p[2]) <= THRESH: edges.append((a,b))
if len(nodes) != N_EXPECT:
    sys.exit(f"!! genome count mismatch: distance file {len(nodes)} vs expected {N_EXPECT}")
parent = {n:n for n in nodes}
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb
for a,b in edges: union(a,b)          # single linkage
cl = collections.defaultdict(list)
for n in nodes: cl[find(n)].append(n)
k = len(cl)
def wilson_upper(n, z=1.959963985):
    d = 1 + z*z/n
    centre = (z*z/(2*n)) / d
    half = (z/d) * math.sqrt((0*(1-0)/n) + z*z/(4*n*n))
    return (centre + half) * 100
with open(OUT_MEMB,"w") as fh:
    fh.write("cluster_id\trepresentative\tn_members\tmembers\n")
    for i,(rep,mem) in enumerate(sorted(cl.items(), key=lambda x:-len(x[1])),1):
        fh.write(f"LP_{i:04d}\t{sorted(mem)[0]}\t{len(mem)}\t{','.join(sorted(mem))}\n")
with open(OUT_SUM,"w") as fh:
    fh.write("species\tn_genomes\tn_clusters\tani_threshold\tmethod\tpositives\tprev_derep\tci_upper_derep\n")
    fh.write(f"Lactiplantibacillus plantarum\t{len(nodes)}\t{k}\t99.5\t"
             f"Mash d<=0.005, single-linkage\t0\t0.00\t{wilson_upper(k):.2f}\n")
print(f"genomes {len(nodes)} -> clusters {k}")
print(f"0/{k}  Wilson upper bound {wilson_upper(k):.2f}%")
print(f"for reference, 0/1665 undereplicated: {wilson_upper(1665):.2f}%")
