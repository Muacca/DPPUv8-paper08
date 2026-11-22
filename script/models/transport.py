"""Numerical massive-shell propagation on the charged Nil history."""
from math import pi, sqrt, cosh, sinh, exp
from itertools import product
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def background(s, charge, holst):
    h = sqrt((1 + holst**2) * charge**2 + 1 / 12)
    d = sqrt(h*h + .25)
    E = exp(2*h*s)
    Y = cosh(2*d*s) - h/d*sinh(2*d*s)
    v = (h*cosh(2*d*s) - d*sinh(2*d*s))/Y
    a2 = E*E*Y
    return E, a2, 1/Y, a2/sqrt(Y), 2*(v-h)


def propagator(field, interval, count):
    # SU(2) evolution of column one; determinant and orthogonality fix column two.
    def rhs(s, flat):
        x, z = field(s)
        p, q = flat.reshape(2, count)
        return -1j*np.r_[z*p+x*q, x*p-z*q]
    y = np.r_[np.ones(count), np.zeros(count)].astype(complex)
    sol = solve_ivp(rhs, interval, y, method="DOP853", rtol=1.5e-12, atol=8e-15)
    assert sol.success
    p, q = sol.y[:, -1].reshape(2, count)
    drift = float(np.max(abs(abs(p)**2+abs(q)**2-1)))
    if drift > 3e-9:
        raise ArithmeticError(drift)
    return p, q, drift


def transitions(p, q, first, last):
    # Project initially negative spinors onto final positive spinors.
    def plus(axis):
        x, y, z = axis
        ang = np.arctan2(np.hypot(x,y), z)
        phase = np.exp(1j*np.arctan2(y, x))
        return np.cos(ang/2), phase*np.sin(ang/2)
    a, b = plus(first)
    c, d = plus(last)
    i, j = -b.conj(), a
    return abs(c*(p*i-q.conj()*j) + d.conj()*(q*i+p.conj()*j))**2


def massive(charge, holst, interval):
    # k = r |j+delta_3|. k<=40 contains the analytic crossing window.
    base = [(f"k{k}",-k,0) for k in range(1,41)]
    base += [("base0",0,0),("base+.5",0,.5),("base-.5",0,-.5)]
    cases=list(product(base, (.25,.5,1.), (1,-1)))
    m = np.array([c[1] for c in cases])
    ks = np.array([c[0][1] for c in cases])
    ls = np.array([c[0][2] for c in cases])
    shift=np.array([1.5*charge*c[2] for c in cases])
    def axes(s, switched=True):
        E,a2,c2,W,_ = background(s,charge,holst)
        return m*W, (-c2/4+ks*a2/(2*pi)+ls*E) + (shift if switched else 0)
    p,q,drift=propagator(axes,interval,len(cases))
    xi,zi=axes(interval[0],False)
    xf,zf=axes(interval[1],False)
    LC=transitions(p,q,(xi,0*xi,zi),(xf,0*xf,zf))
    LR=abs(q)**2
    records=[]
    for i in range(0,len(cases),2):
        label=cases[i][0][0]
        k=-ks[i]
        if k>6:
            # Keep all possibly crossing levels, plus three of each fiber type.
            h=sqrt((1+holst**2)*charge**2+1/12)
            d=sqrt(h*h+.25)
            if k > (1.5*charge) * (2*pi) * (2*d) / exp(4*h*interval[0]):
                continue
        roots = {}
        grid=np.linspace(*interval,81)
        for j,sign in enumerate(("plus","minus")):
            def f(s):
                return axes(s)[1][i+j]
            roots[sign]=[brentq(f,a,b) for a,b in zip(grid,grid[1:]) if f(a)*f(b)<0]
        def measure(v):
            a,b=v[i:i+2]
            return dict(plus=a,minus=b,diff=a-b,ratio=a/b if b else None)
        records.append(dict(label=label,mass=m[i],LC=measure(LC),LR=measure(LR),roots=roots))
    return dict(u=charge,h=holst,I=interval,unitarity=drift,rows=records)


