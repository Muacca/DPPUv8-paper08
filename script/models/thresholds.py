"""Full-line signed Hermite residuals plus min-max isolation at the four target periods."""
import json
from math import ceil,floor

from flint import arb,ctx
import numpy as np
from scipy.linalg import eigh

from . import hermite
ctx.prec=170


def line(tr,n):
    disc=arb(tr*tr-4).sqrt()
    b=2*arb.pi()*(2*n/disc).sqrt()
    # The original function's b=2*pi*sqrt(argument)/3^(1/4).
    argument=2*arb(3).sqrt()*n/disc
    count=20
    D,H=hermite.forms(argument,count)
    trial=eigh(np.array([[float(e) for e in row] for row in H]),subset_by_index=[0,1])[1].T
    spectrum=[]
    for vec in trial:
        v=[round(z*10**9) for z in vec]
        nonzero=[(i,x) for i,x in enumerate(v) if x]
        norm=sum(x*x for x in v)
        def form(M):
            return sum((x*y*M[i][j] for i,x in nonzero for j,y in nonzero),arb(0))/norm
        mean=form(D);residual=(form(H)-mean**2).sqrt()
        ev=mean+arb(0,residual.upper())
        spectrum.append((ev,dict(mean=str(mean),residual=str(residual),integer_trial=v)))
    delta=arb(7)/50
    omega=b*(2*(1-delta)).sqrt()
    second=b*b+b+omega-1/(4*delta)
    third=b*b-b+3*omega-1/(4*delta)
    first=spectrum[0][0];exc=spectrum[1][0]
    assert first>0 and exc<0 and first**2<second and exc**2<third and first**2<exc**2
    certificate={"n_abs_F":n,"b":str(b),"second_squared_lower":str(second),"third_squared_lower":str(third),
                 "hermite_count":count,"width":7,"positive_product_states":[z[1] for z in spectrum]}
    return first,-exc,certificate


def calculate(p):
    ctx.prec=170
    rows=json.loads((p/"classes.json").read_bytes())["classes"]
    output=[]
    for row in rows:
        if row["period"] not in ([3],[2,3],[2,3,4],[2,2,4,4]):
            continue
        tr=row["trace"]; N=row["norm_F"]
        ell=((tr+arb(tr*tr-4).sqrt())/2).log()
        sectors={-1:(N["positive"],len(N["positive_reps"])),1:(N["negative_abs"],len(N["negative_reps"]))}
        lines={n:line(tr,n) for n,mult in sectors.values()}
        records={}
        for sign,(n,mult) in sectors.items():
            other,other_mult=sectors[-sign]
            ground=lines[n][0];candidate=lines[other][1]
            if ground<candidate:
                val=ground;origin="ground";dim=mult;source=n
            else:
                assert candidate<ground
                val=candidate;origin="first_excited";dim=other_mult;source=other
            charge=4*ell*val
            # Strict outward decimal endpoints verified again in ball arithmetic.
            left=floor(float(charge.lower())*1e5)-1
            right=ceil(float(charge.upper())*1e5)+1
            assert arb(left)/100000<charge and charge<arb(right)/100000
            records["+" if sign==1 else "-"]={"abs_Q_enclosure":[left/100000,right/100000],"origin":origin,
                "complex_dimension":dim,"n_abs_F":source,"other_candidate_abs_Q":str(4*ell*(candidate if origin=="ground" else ground))}
        output.append({"period":row["period"],"trace":tr,"levels":[z[2] for z in lines.values()],"results":records})
    (p/"thresholds.json").write_text(json.dumps(output,indent=2)+"\n",encoding="ascii")
    for row in output:
        print("period",row["period"],json.dumps(row["results"]))


