"""Character multiplicities and ordinary eta on fifteen spin lifts."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from flint import arb
import math

# (signed class count, SU2 eigenangle theta/(2*pi)).
tet=[(1,F(0)),(1,F(1,2)),(6,F(1,4)),(8,F(1,6)),(8,F(1,3))]
oct_extra=[(12,F(1,4)),(6,F(1,8)),(6,F(3,8))]
ico=[(1,F(0)),(1,F(1,2)),(30,F(1,4)),(20,F(1,6)),(20,F(1,3))]+[(12,F(j,10)) for j in range(1,5)]
def chars(n,alpha):
    if alpha==0:return arb(n+1)
    if alpha==F(1,2):return arb((-1)**n*(n+1))
    t=2*arb.pi()*alpha.numerator/alpha.denominator
    return (t*((n+1)%alpha.denominator)).sin()/t.sin()
def binary(order,classes,right=False):
    @lru_cache(None)
    def invariant(n):
        ball=sum((w*chars(n,a) for w,a in classes),arb(0))/order
        i=ball.unique_fmpz()
        assert i is not None and ball.rad()<arb('1e-30')
        return int(i)
    def shell(k):
        p=(k+2)*invariant(k);m=(k+1)*invariant(k+1)
        return (m,p) if right else (p,m)
    return shell
def lens(p,a,b):
    def count(l,r):
        h=Counter(a*u%p for u in range(-l,l+1,2))
        return sum(h[-b*v%p] for v in range(-r,r+1,2))
    return lambda k:(count(k,k+1),count(k+1,k))
def cases():
    return [
        ('S3',1,lens(1,0,0)),('RP3_even',2,lens(2,1,0)),('RP3_odd',2,lens(2,0,1)),
        ('C3_left',3,lens(3,1,0)),('C3_right',3,lens(3,0,1)),
        ('L5',5,lens(5,2,1)),('L7_23',7,lens(7,3,2)),
        ('2T',12,binary(24,tet)),('2O',24,binary(48,tet+oct_extra)),
        ('2O_sign',24,binary(48,tet+[(-w,a) for w,a in oct_extra])),('2I',60,binary(120,ico)),
        ('Q8',4,binary(8,[(1,F(0)),(1,F(1,2)),(6,F(1,4))]))]+[
        ('Q8_'+axis,4,binary(8,[(1,F(0)),(1,F(1,2)),(-2,F(1,4))])) for axis in ('i','j','k')]
def coefficients(P,fn):
    d=[p-m for p,m in [fn(k) for k in range(4*P)]]
    polys=[]
    for r in range(P):
        a,b,c,e=d[r::P]
        A=F(c-2*b+a,2*P*P)
        B=F(b-a,P)-A*(2*r+P)
        C=a-B*r-A*r*r
        assert e==A*(r+3*P)**2+B*(r+3*P)+C
        polys.append([C,B,A])
    assert not any(sum(col) for col in zip(*polys))
    num=[sum((-1)**j*math.comb(3,j)*d[i-j*P] for j in range(i//P+1)) for i in range(3*P)]
    for m in range(3):
        assert sum(math.comb(i,m)*v for i,v in enumerate(num) if i>=m)==0
    eta=F(-sum(math.comb(i,3)*v for i,v in enumerate(num) if i>=3),P**3)
    return eta,polys,max((k for row in polys for k,x in enumerate(row) if x),default=0)
