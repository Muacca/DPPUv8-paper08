"""Exact lowest shells in units of 2*pi, using the horizontal shortest vectors."""
from fractions import Fraction as F
from itertools import product


orders = [1,2,3,4,6,4]
sizes = ["1/2","1","2"]

def spins(g):
    if g == 1:
        yield from product((0,1), repeat=3)
    elif g == 6:
        yield from product((1,-1), repeat=2)
    else:
        horizontal = {2:list(product((0,1), repeat=2)), 3:[(0,0)],
                      4:[(0,0),(1,1)], 5:[(0,0)]}[g]
        yield from ((s,i,j) for s in [1,-1] for i,j in horizontal)

def horizontal_min(i,j,hexagon):
    if hexagon:
        assert (i,j) == (0,0)
        return F(4,3),6
    return (F(i+j,4),2**(i+j)) if i+j else (F(1),4)

def result(g, spin):
    m = orders[g-1]
    if g == 1:
        i,j,k = spin
        axial = i == j == 0
    elif g == 6:
        i=j=k=1
        axial=False
    else:
        sign,i,j = spin
        k = int(-sign**m < 0)
        axial = i==j==0
    r,n = horizontal_min(i,j,g in (3,5))
    z = F(k,2)
    eta = F(0)
    zero = 2 if g==1 and spin==(0,0,0) else 0
    crit = {}
    if axial and g not in (1,6):
        c = (F(-1,2) + (F(m,2) if sign<0 else 0)) % m
        eta = 2*(1-2*c/m)
        crit = dict(split_L_squared=str((min(c,m-c)**2-z*z)/r),
                    both_axial_L_squared=str((max(c,m-c)**2-z*z)/r))
    rows = {}
    for label in sizes:
        L = F(label)
        transverse = r+(z/L)**2
        copies = F(n*(2 if z else 1),m)
        assert copies.denominator==1
        pair=[]
        for branch in (1,-1):
            choices=[(transverse,int(copies),"H")]
            if axial:
                height = (z or 1) if g==1 else (c if branch==1 else m-c)
                choices.append(((height/L)**2,2,"A"))
            level = min(row[0] for row in choices)
            winners = [row for row in choices if row[0]==level]
            pair.append(dict(lambda_squared_in_2pi_units=str(level),
                             dim=sum(v[1] for v in winners),origin="+".join(v[2] for v in winners)))
        rows[label] = pair
    return dict(group=g,spin=list(spin),eta=str(eta),lc_kernel=zero,critical=crit,thresholds=rows)

def calculate():
    out=[result(g,spin) for g in range(1,7) for spin in spins(g)]
    assert len(out)==28 and sum(F(row["eta"])!=0 for row in out)==8
    assert [sum(F(v["eta"])!=0 for v in out if v["group"]==g) for g in range(1,7)]==[0,2,2,2,2,0]
    assert sum(row["lc_kernel"] for row in out)==2

    return out
