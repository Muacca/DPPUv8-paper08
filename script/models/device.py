"""Validated midpoint SU(2) propagation for the two representative charge pairs."""
from flint import arb, acb, ctx


ctx.prec=90


def certificate(sign,holst):
    ctx.prec=90
    u=arb(3)/8
    H=((1+holst**2)*u**2+arb(1)/12).sqrt()
    d=(H**2+arb(1)/4).sqrt()
    nu=1/(2*arb.pi())
    def field(t):
        E=(2*H*t).exp()
        Y=(2*d*t).cosh()-H/d*(2*d*t).sinh()
        v=(H*(2*d*t).cosh()-d*(2*d*t).sinh())/Y
        vp=-1/(2*Y**2)
        A=E**2*Y
        x=E**2*Y.sqrt()/4
        g=1/(4*Y)
        z=-g-nu*A+sign*3*u/2
        xp=x*(4*H-v)
        zp=-2*v*g-nu*A*(4*H-2*v)
        xpp=x*((4*H-v)**2-vp)
        zpp=-g*(4*v**2+2*vp)-nu*A*((4*H-2*v)**2-2*vp)
        return x,z,(abs(xpp)+abs(zpp))/24+(abs(x)+abs(z))*(abs(xp)+abs(zp))/6
    def positive_lc(t):
        x,z,_=field(t)
        z-=sign*3*u/2
        norm=(x**2+z**2).sqrt()
        return ((1+z/norm)/2).sqrt(),((1-z/norm)/2).sqrt()
    steps=4096
    step=arb(1.5)/steps
    bound=arb(0)
    p,q=acb(1),acb(0)
    I=acb(0,1)
    for j in range(steps):
        mid=(j+arb(.5))*step
        x,z,_=field(mid)
        om=(x**2+z**2).sqrt()
        C=(om*step).cos()
        S=-I*(om*step).sin()/om
        p,q=(C+S*z)*p+S*x*q,S*x*p+(C-S*z)*q
        _,_,local=field(arb(mid,step/2))
        bound+=local.upper()*step**3
    i,j=positive_lc(arb(0))
    c,e=positive_lc(arb(1.5))
    amp=c*(-j*p-i*q.conjugate())+e*(-j*q+i*p.conjugate())
    error=2*bound+bound**2
    return dict(LC=arb(abs(amp)**2,error.upper()),LR=arb(abs(q)**2,error.upper()))


