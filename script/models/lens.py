"""Exact Abel limit from weight counts, with independent Arb group averages."""
from fractions import Fraction as F
from math import gcd, comb
from collections import Counter
from flint import arb, ctx


ctx.prec = 160

def rank(p,a,b,l,r):
    hist = Counter(a*i%p for i in range(-l,l+1,2))
    return sum(hist[-b*j%p] for j in range(-r,r+1,2))

def classify(p,r,t):
    if p<2 or gcd(p,r)!=1 or t not in (0,1) or (1+r+t*p)%2:
        raise ValueError("invalid space or lift")
    a,b=(1+r+t*p)//2,(1-r+t*p)//2
    d=[rank(p,a,b,k,k+1)-rank(p,a,b,k+1,k) for k in range(3*p)]
    # (1-z^p)^3 times the signed generating function.
    coeff = [sum((-1)**j*comb(3,j)*d[k-j*p] for j in range(k//p+1)) for k in range(3*p)]
    assert all(sum(comb(k,s)*v for k,v in enumerate(coeff) if k>=s)==0 for s in range(3))
    eta=F(-sum(comb(k,3)*v for k,v in enumerate(coeff) if k>=3),p**3)
    mp,mm=rank(p,a,b,0,1),rank(p,a,b,1,0)
    assert (mp,mm)==(2*int(b%p==0),2*int(a%p==0))
    left,right=3*eta,5*eta-2*d[0]
    z=3*d[0]/(2*(d[0]-eta)) if d[0]!=eta else None
    kind='identical' if left==right==0 else 'isolated' if left*right<0 else 'none'
    z=z if kind=='isolated' else ''
    return dict(p=p,r=r,t=t,a=a%p,b=b%p,eta=eta,m_plus=mp,m_minus=mm,
                open_type=kind,open_zero=z,left=left,right=right,
                left_zero=left==0,right_zero=right==0)

def group_check(row):
    p,a,b,r=row['p'],row['a'],row['b'],row['r']
    eta=arb(0)
    pos=arb(2)
    neg=arb(2)
    for k in range(1,p):
        theta=2*arb.pi()*k/p
        l,h=2*(a*theta).cos(),2*(b*theta).cos()
        pos+=h
        neg+=l
        eta+=2*(h-l)/((2-2*theta.cos())*(2-2*(r*theta).cos()))
    for value,target in [(pos/p,F(row['m_plus'])),(neg/p,F(row['m_minus'])),(eta/p,row['eta'])]:
        assert value.contains(arb(target.numerator)/target.denominator) and value.rad()<arb('1e-30')

