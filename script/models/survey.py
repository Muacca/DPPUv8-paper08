"""Enumerate primitive modular necklaces and exact covector-norm minima."""
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt
import json


UNIT=(1,0,0,1)
L=(1,1,0,1)
R=(1,0,1,1)


def mul(X,Y):
    a,b,c,d=X
    e,f,g,h=Y
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


def canon(seq):
    return min(seq[j:]+seq[:j] for j in range(len(seq)))


def product(factors):
    A=UNIT
    for M in factors:
        A=mul(A,M)
    return A


def inverse(A):
    a,b,c,d=A
    det=a*d-b*c
    assert abs(det)==1
    return det*d,-det*b,-det*c,det*a


def symmetry_witness(word,reverse):
    swap=word.translate(str.maketrans("LR","RL"))
    transformed=swap[::-1] if reverse else swap
    transform=(0,-1,1,0) if reverse else (0,1,1,0)
    for k in range(len(word)):
        if transformed[k:]+transformed[:k]==word:
            prefix=product(L if c=="L" else R for c in transformed[:k])
            # Return in the integral basis of the negative continued fraction.
            return mul(L,mul(inverse(prefix),mul(transform,inverse(L))))
    return None


def necklaces(limit):
    # Each mixed circular word can start LR. Every further factor strictly increases trace.
    stack=[("LR",mul(L,R))]
    found=set()
    while stack:
        word,A=stack.pop()
        if A[0]+A[3]>limit:
            continue
        if word not in (word+word)[1:-1]:
            found.add(canon(word))
        stack.extend((word+s,mul(A,m)) for s,m in (("L",L),("R",R)))
    return sorted(found)


def surdsign(a,b,D):
    if a*b>=0:
        return (a+b>0)-(a+b<0)
    return ((a>0)-(a<0))*((a*a>D*b*b)-(a*a<D*b*b))


def orbit_minima(M):
    a,b,c,d=M
    trace=a+d;disc=trace**2-4
    assert b>0 and c>0
    maximum=max(b,c)
    ymax=isqrt(4*b*maximum*trace//disc)
    smax=isqrt(4*b*maximum*trace)
    bins={}
    for y in range(-ymax,ymax+1):
        offset=(d-a)*y
        xmin=-((smax+offset)//(2*b))
        xmax=(smax-offset)//(2*b)
        for x in range(xmin,xmax+1):
            n=b*x*x+(d-a)*x*y-c*y*y
            if not 0<abs(n)<=maximum:
                continue
            # lambda^-1 <= (s+sqrt(D)y)^2/(4 b |n|) < lambda.
            s=2*b*x+offset
            base=2*(s*s+disc*y*y)-4*b*abs(n)*trace
            cross=4*s*y
            if surdsign(base,cross+4*b*abs(n),disc)>=0 and surdsign(base,cross-4*b*abs(n),disc)<0:
                bins.setdefault(n,[]).append([x,y])
    pos=min(n for n in bins if n>0)
    neg=min(-n for n in bins if n<0)
    return {"positive":pos,"negative_abs":neg,"positive_reps":bins[pos],"negative_reps":bins[-neg],
            "search_cap":maximum,"smax":smax,"ymax":ymax}


def single(word):
    # Rotate so the last letter is R, then each L-run and its next R gives one digit.
    j=word.index("R")+1
    shifted=word[j:]+word[:j]
    digits=canon(tuple(2+len(run) for run in shifted[:-1].split("R")))
    aligned="".join("L"*(v-2)+"R" for v in digits)
    positive=product(L if s=="L" else R for s in aligned)
    A=product((v,-1,1,0) for v in digits)
    assert mul(A,L)==mul(L,positive)
    tr=A[0]+A[3]
    comm=symmetry_witness(aligned,False)
    rev=symmetry_witness(aligned,True)
    U=comm is not None
    reversal=rev is not None
    if U:
        assert mul(comm,A)==mul(A,comm) and comm[0]*comm[3]-comm[1]*comm[2]==-1
    if reversal:
        assert mul(rev,A)==mul(inverse(A),rev) and rev[0]*rev[3]-rev[1]*rev[2]==1
    norm=orbit_minima(positive)
    sigma=sum(v-3 for v in digits)
    return {"trace":tr,"period":digits,"A":A,"positive_word":aligned,"positive_A":positive,
            "type_U":U,"type_R":reversal,"chiral":not (U or reversal),"eta_D":str(Fraction(sigma,6)),
            "U_matrix":comm,"R_matrix":rev,
            "form_content":gcd(positive[1],gcd(positive[2],positive[0]-positive[3])),
            "norm_F":norm,"thresholds_split":norm["positive"]!=norm["negative_abs"]}


def calculate(out):
    rows=sorted((single(w) for w in necklaces(60)),key=lambda row:(row["trace"],row["period"]))
    summary=Counter()
    for row in rows:
        summary["classes"]+=1
        summary["chiral"]+=row["chiral"]
        summary["type_R_only"]+=row["type_R"] and not row["type_U"]
        summary["type_U_only"]+=row["type_U"] and not row["type_R"]
        summary["type_both"]+=row["type_R"] and row["type_U"]
        summary["chiral_eta_zero"]+=row["chiral"] and row["eta_D"]=="0"
        summary["split"]+=row["thresholds_split"]
    (out/"classes.json").write_text(json.dumps({"summary":summary,"classes":rows},separators=(",",":"))+"\n",encoding="ascii")
    table=["# Classes with trace 3 through 60","",
           "The signed minima are the magnitudes of the integral form F. Physical km has the opposite sign.","",
           "| trace | period | U | R | Dirac eta | min F>0 | min -F>0 | split |",
           "|---:|---|:---:|:---:|---:|---:|---:|:---:|"]
    for r in rows:
        n=r["norm_F"]
        signs=["yes" if r[k] else "no" for k in ("type_U","type_R","thresholds_split")]
        table.append(f'| {r["trace"]} | {",".join(map(str,r["period"]))} | {signs[0]} | {signs[1]} | {r["eta_D"]} | {n["positive"]} | {n["negative_abs"]} | {signs[2]} |')
    (out/"all_classes.md").write_text("\n".join(table)+"\n",encoding="ascii")
    print(json.dumps(summary,sort_keys=True))
