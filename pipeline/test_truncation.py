import os,subprocess,collections,random
random.seed(42)
REF="refs/gate2_reference.faa"; OUT="search/trunc_test"
os.makedirs(OUT,exist_ok=True)
seqs={}; lab=None; buf=[]
for line in open(REF):
    if line[0]==">":
        if lab: seqs[lab]="".join(buf)
        lab=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
if lab: seqs[lab]="".join(buf)
tgt=[k for k in seqs if k.startswith("TARGET_tyrDC")]
if not tgt: raise SystemExit("TARGET_tyrDC 참조가 없습니다")
base=tgt[0]; S=seqs[base]; L=len(S)
print("기준 서열:",base,L,"aa")
AA="ACDEFGHIKLMNPQRSTVWY"
cases=[]
for frac in (0.95,0.90,0.85,0.80,0.70,0.50):
    cases.append(("TRUNC_%d"%(frac*100), S[:int(L*frac)]))
for pct in (5,10,20,40):
    s=list(S)
    for i in random.sample(range(L), int(L*pct/100)): s[i]=random.choice(AA)
    cases.append(("POINT_%d"%pct, "".join(s)))
cases.append(("INTACT", S))
q=os.path.join(OUT,"query.faa")
with open(q,"w") as f:
    for n,s in cases:
        f.write(">SYNTH_%s|%s synthetic control\n"%(n,n))
        for i in range(0,len(s),60): f.write(s[i:i+60]+"\n")
print("합성 대조",len(cases),"건")
h=os.path.join(OUT,"hits.tsv")
subprocess.run(["diamond","blastp","-q",q,"-d","search/ref","-o",h,"--outfmt","6",
  "qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore",
  "--evalue","1e-5","--max-target-seqs","20","--quiet"],check=True)
ID_P,COV_P,MAR=60.0,80.0,10.0; ID_A,COV_A=40.0,50.0; TR=0.80
per=collections.defaultdict(list)
for line in open(h):
    c=line.rstrip("\n").split("\t")
    if len(c)<8: continue
    per[c[0]].append({"fam":c[1].split("|")[0],"pident":float(c[2]),"length":int(c[3]),
                      "qlen":int(c[4]),"slen":int(c[5]),"bits":float(c[7])})
print("")
print("%-12s %7s %7s %7s %8s  %s"%("case","pident","cov%","qlen/sl","margin","call / variant"))
for n,_ in cases:
    key="SYNTH_%s|%s"%(n,n)
    hits=per.get(key,[])
    if not hits:
        print("%-12s  히트 없음"%n); continue
    byf=collections.defaultdict(list)
    for x in hits: byf[x["fam"]].append(x)
    rk=sorted(((f,max(v,key=lambda y:y["bits"])) for f,v in byf.items()),key=lambda z:-z[1]["bits"])
    fam,top=rk[0]; second=rk[1][1]["bits"] if len(rk)>1 else 0.0
    mg=top["bits"]-second
    cov=100.0*top["length"]/max(top["slen"],1); lr=top["qlen"]/max(top["slen"],1)
    if not fam.startswith("TARGET"): call,var="absent","decoy_top"
    elif top["pident"]>=ID_P and cov>=COV_P and mg>=MAR:
        var="truncated" if lr<TR else ("point-variant" if top["pident"]<95.0 else "intact")
        call="present"
    elif top["pident"]>=ID_A and cov>=COV_A: call,var="ambiguous",""
    else: call,var="absent",""
    print("%-12s %7.1f %7.1f %7.2f %8.1f  %s / %s"%(n,top["pident"],cov,lr,mg,call,var))
print("")
print("truncated 가 한 번도 안 나오면 그 분류는 발동 불가능하다는 뜻입니다.")
