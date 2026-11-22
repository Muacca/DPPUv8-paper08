"""Recompute Paper 08 tables and certificates, or render its comparison figure."""
import argparse
import csv
import json
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=True, allow_nan=False, indent=1) + "\n", encoding="utf-8")


def recompute(out):
    from flint import ctx, arb
    from models import survey, thresholds, spheres, lens, nil, flat, device, transport

    survey.calculate(out)
    thresholds.calculate(out)
    ctx.prec = 160
    rows = []
    for name, period, shell in spheres.cases():
        eta, polys, degree = spheres.coefficients(period, shell)
        mp, mm = shell(0)
        rows.append(dict(name=name, period=period, degree=degree, eta=str(eta),
                         m_plus=mp, m_minus=mm, ascending_polynomials=[[str(c) for c in p] for p in polys]))
    save(out / "shells.json", rows)
    cases = []
    for p in range(2, 33):
        for r in range(1, p):
            if gcd(r, p) != 1:
                continue
            for t in (0, 1):
                if (1 + r + t*p) % 2:
                    continue
                row = lens.classify(p, r, t)
                lens.group_check(row)
                cases.append(row)
    assert len(cases) == 434
    with (out / "lens.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    ctx.prec = 192
    leaves = []
    for r in range(1, 5):
        for spin in product((0, 1), (0, 1), range(1 + int(r % 2 == 0))):
            for t in (0.5, 1, 2):
                eta, jumps, kernel = nil.eta(r, spin, arb(1), arb(t), arb(0))
                record = nil.thresholds(r, spin, arb(t))
                assert kernel == record["full"]["zero_kernel"]
                leaves.append(dict(r=r, spin2=spin, t=t, eta=nil.ball(eta), J=jumps, **record))
    save(out / "nil.json", leaves)
    save(out / "flat.json", {"spins": flat.calculate()})
    numerical = [transport.massive(u, h, interval)
                 for u, h, interval in product((.125, .375), (0, 1), ((-.5, 0.), (0., 1.5)))]
    save(out / "nil_measurements.json", numerical)
    certificates = []
    for h in (0, 1):
        pos, neg = device.certificate(1, h), device.certificate(-1, h)
        for name in ("LC", "LR"):
            assert pos[name] > neg[name]
            certificates.append(dict(h=h, device=name, plus=str(pos[name]), minus=str(neg[name]),
                                     difference=str(pos[name]-neg[name]), ratio=str(pos[name]/neg[name])))
    save(out / "device_certificate.json", certificates)
    return verify(out)


def read(folder, name):
    return json.loads((folder / name).read_text(encoding="utf-8-sig"))


def verify(folder):
    import math
    from flint import arb
    classes = read(folder, "classes.json")
    c = classes["classes"]
    assert len(c) == len({r["positive_word"] for r in c}) == 488
    assert sum(r["chiral"] for r in c) == sum(r["thresholds_split"] for r in c) == 472
    assert all(r["chiral"] == r["thresholds_split"] for r in c)
    for r in c:
        a,b,cc,d = r["positive_A"]
        assert a*d-b*cc==1 and a+d==r["trace"]
        for side, sign in (("positive", 1), ("negative", -1)):
            target = sign*r["norm_F"]["positive" if sign==1 else "negative_abs"]
            for x,y in r["norm_F"][side+"_reps"]:
                assert b*x*x+(d-a)*x*y-cc*y*y == target
        assert Fraction(r["eta_D"]) == Fraction(r["positive_word"].count("L")-r["positive_word"].count("R"), 6)
    assert sum(r["chiral"] and Fraction(r["eta_D"])==0 for r in c) == 6
    assert sum(r["type_R"] and not r["type_U"] for r in c)==8
    golden_period=math.log((3+math.sqrt(5))/2)
    assert 22.9843 < 4*golden_period*math.sqrt(35.6460648209341)
    assert 4*golden_period*math.sqrt(37.8505887358739) < 23.6845
    beta=1.5*math.log(math.tan(math.pi/6)/math.tan(math.pi/12))
    assert math.isclose(math.cos(beta),.4071382691,abs_tol=5e-10)
    assert math.isclose(math.sin(beta),.9133665364,abs_tol=5e-10)
    spherical=read(folder, "shells.json")
    assert len(spherical)==15
    assert {r["name"]:r["eta"] for r in spherical}["L7_23"]=="4/7"
    with (folder / "lens.csv").open(encoding="utf-8-sig",newline="") as handle:
        lenses=list(csv.DictReader(handle))
    assert len(lenses)==434 and {r["open_zero"] for r in lenses if r["open_type"]=="isolated"}=={"12/7","27/14","24/11"}
    flat=read(folder,"flat.json")["spins"]
    assert len(flat)==28 and sum(bool(Fraction(r["eta"])) for r in flat)==8
    leaves=read(folder,"nil.json")
    assert len(leaves)==72
    assert all(not arb(r["nonzero"]["Q_positive"]["abs_Q"]).overlaps(arb(r["nonzero"]["Q_negative"]["abs_Q"])) for r in leaves if r["t"]==1)
    sol=read(folder,"thresholds.json")
    assert len(sol)==4 and next(r for r in sol if r["period"]==[2,3])["results"]["+"]["origin"]=="first_excited"
    measurements=read(folder,"nil_measurements.json")
    probs=next(r for r in measurements if r["u"]==.375 and r["h"]==0 and list(r["I"])==[0,1.5])
    obs=next(r for r in probs["rows"] if r["label"]=="k1" and r["mass"]==.25)["LC"]
    cert=next(r for r in read(folder,"device_certificate.json") if r["h"]==0 and r["device"]=="LC")
    assert all(arb(cert[k]).contains(arb(obs[k])) for k in ("plus","minus","ratio"))
    assert math.isclose(obs["plus"],.168451,abs_tol=5e-7) and math.isclose(obs["minus"],.042475,abs_tol=5e-7)
    print("Verified 488 Sol classes, four whole-line certificates, 15 spherical lifts, 434 lens lifts, 28 flat lifts, 72 Nil leaves and the two-charge detector.")
    return True


def render(folder, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = read(folder, "thresholds.json")
    n = next(r for r in read(folder, "nil_measurements.json") if r["u"] == .375 and r["h"] == 0 and r["I"] == [0, 1.5])
    rows = [r for r in n["rows"] if r["label"]=="k1"]
    fig, (ax,bx) = plt.subplots(1,2,figsize=(7.3,3.15), layout="constrained")
    colors = ("#11486C","#A54523")
    for j,(sgn,label) in enumerate((("-",r"$Q_\Sigma<0$"),("+",r"$Q_\Sigma>0$"))):
        for i,row in enumerate(t):
            lo,hi=row["results"][sgn]["abs_Q_enclosure"]
            ax.plot((lo+hi)/2,i+(j-.5)*.17,"os"[j],color=colors[j],label=label if i==0 else "")
        bx.plot([r["mass"] for r in rows], [r["LC"]["plus" if j==0 else "minus"] for r in rows],
                marker="os"[j], color=colors[j], label=label)
    ax.set(yticks=range(len(t)), yticklabels=[str(tuple(r["period"])).replace(",)",")") for r in t], xlabel="First charge magnitude", ylabel="Sol period")
    bx.set(xlabel="Dirac mass",ylabel="Mean occupation per copy",xticks=[.25,.5,1],ylim=(0,.19))
    ax.set_title("Nonzero fibre sector",fontsize=10)
    bx.set_title("Same endpoint detector on Nil",fontsize=10)
    for a in (ax,bx):
        a.grid(axis="x",alpha=.18)
        a.spines[["top","right"]].set_visible(False)
        a.legend(fontsize=8,frameon=False)
    fig.savefig(out / "readouts.pdf",metadata={"CreationDate":None,"ModDate":None})
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="mode",required=True)
    for name in ("verify","recompute","render"):
        cmd=sub.add_parser(name)
        if name!="recompute":
            cmd.add_argument("--data",type=Path,default=Path(__file__).resolve().parents[1]/"data")
        if name!="verify":
            cmd.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.mode=="verify":
        verify(args.data)
        return
    args.output.mkdir(parents=True,exist_ok=False)
    if args.mode=="recompute":
        recompute(args.output)
    else:
        render(args.data,args.output)


if __name__=="__main__":
    main()
