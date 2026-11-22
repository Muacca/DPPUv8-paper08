"""Heisenberg spatial spectra and fixed-leaf charge thresholds."""
from collections import Counter
from flint import arb, ctx

ctx.prec = 192
PI = arb.pi()
def sign(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    assert x == 0, f"unresolved sign: {x}"
    return 0


def upper(x):
    # Ceiling of a strict certified upper bound; an extra integer is harmless.
    return int(x.upper().ceil().unique_fmpz()) + 1


def ball(x):
    return x.str(21)


def base_shells(x, y, bound):
    return Counter((2*i+x)**2 + (2*j+y)**2
                   for i in range(-bound, bound+1)
                   for j in range(-bound, bound+1))


def eta(r, spin, a, c, u):
    x, y, z = spin
    t = c/(a*a)
    vol0 = (2*PI)**4 / r
    w = a*a*c
    mu = 3*u/(2*w)
    offset = mu - t/4
    direction = sign(offset)
    gap = abs(offset)
    hurwitz = arb(r)/6 if z == 0 else -arb(r)/12
    smooth = vol0*w/PI**2 * (((t/4)**3-mu**3)/3-offset*t**2/48)
    jump = kernel = 0
    if z == 0:
        for norm4, mult in base_shells(x, y, upper(a*gap)+1).items():
            if norm4 == 0:
                jump += 2*direction
                kernel += 2*(direction == 0)
            else:
                cmp = sign(gap - arb(norm4).sqrt()/(2*a))
                jump += mult*direction*(cmp+1)
                kernel += mult*(cmp == 0)
    for k in range(r if z == 0 else r//2, upper(2*PI*c*gap)+1, r):
        nu = k/(2*PI)
        if direction > 0:
            cmp = sign(gap - nu/c)
            jump += 2*k*(cmp+1)
            kernel += 2*k*(cmp == 0)
        # Pair radius grows strictly with n. The vertical term alone bounds k.
        for n in range(1, upper((a*gap)**2/(2*nu))+1):
            cmp = sign(gap**2 - (nu/c)**2 - 2*nu*n/a**2)
            if cmp < 0:
                break
            jump += 2*k*direction*(cmp+1)
            kernel += 2*k*(cmp == 0)
    return hurwitz + smooth + jump, jump, kernel


def thresholds(r, spin, c):
    x, y, z = spin
    k = r if z == 0 else r//2
    nu = k/(2*PI)
    radius = ((nu/c)**2 + 2*nu).sqrt()
    assert radius > c/4
    # Work with eigenvalues at zero shift: Q = 4 V W lambda_LC.
    positive = (radius-c/4, f"P{k}", 2*k)
    negative = (nu/c+c/4, f"U{k}", 2*k)
    out = {}
    for sector in ("nonzero", "full"):
        kernel = 0
        if sector == "full" and z == 0:
            bound = upper(max(positive[0], negative[0])+c/4)+1
            for norm4, mult in base_shells(x, y, bound).items():
                for side in (-1, 1):
                    value = side*arb(norm4).sqrt()/2-c/4
                    sgn = sign(value)
                    if not sgn:
                        kernel += mult
                        continue
                    best = positive if sgn > 0 else negative
                    dist = abs(value)
                    if dist < best[0]:
                        best = (dist, f"B{norm4}{'+' if side>0 else '-'}", mult)
                    elif not dist > best[0]:
                        assert dist == best[0]
                        best = (dist, best[1], best[2]+mult)
                    if sgn > 0:
                        positive = best
                    else:
                        negative = best
        record = {"zero_kernel": kernel}
        for name, entry in (("Q_positive", positive), ("Q_negative", negative)):
            lam, mode, mult = entry
            density = 2*c*lam/3
            record[name] = dict(abs_u=ball(density), abs_Q=ball(6*(2*PI)**4/r*density),
                                dim=mult, mode=mode)
        out[sector] = record
    return out


