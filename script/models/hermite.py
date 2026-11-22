"""Full-line Hermite forms at Gaussian width seven."""
from math import factorial
from flint import arb

def moments(i,j,t,w):
    z=arb(t)/(2*w).sqrt()
    v=sum((z**(i+j-2*l)/(factorial(i-l)*factorial(j-l)*factorial(l))
           for l in range(min(i,j)+1)),arb(0))
    return (arb(t*t)/(4*w)).exp()*arb(factorial(i)*factorial(j)).sqrt()*v


def forms(n, count):
    # In a fixed spin basis C=sigma3 and (-sigma2_old+sigma3_old)/sqrt(2)=sigma1.
    # D = -i sigma2 d/dx + b cosh(x) sigma3 + b sinh(x) sigma1 for positive product.
    w=arb(7)
    b=2*arb.pi()*arb(n).sqrt()/arb(3).sqrt().sqrt()
    size=2*count
    D=[[arb(0) for _ in range(size)] for _ in range(size)]
    H=[[arb(0) for _ in range(size)] for _ in range(size)]
    for i in range(count):
        for j in range(count):
            e=moments(i,j,1,w)
            ch=e if (i+j)%2==0 else arb(0)
            sh=e if (i+j)%2 else arb(0)
            ch2=moments(i,j,2,w) if (i+j)%2==0 else arb(0)
            der=(w*j/2).sqrt()*int(i==j-1)-(w*(j+1)/2).sqrt()*int(i==j+1)
            lap=(w*(2*j+1)/2*int(i==j)-w*arb(j*(j-1)).sqrt()/2*int(i==j-2)
                 -w*arb((j+1)*(j+2)).sqrt()/2*int(i==j+2))
            for spin in (0,1):
                s=1-2*spin
                D[i+spin*count][j+spin*count]=s*b*ch
                D[i+spin*count][j+(1-spin)*count]=b*sh-s*der
                H[i+spin*count][j+spin*count]=lap+b*b*ch2-s*b*ch
                H[i+spin*count][j+(1-spin)*count]=b*sh
    return D,H


