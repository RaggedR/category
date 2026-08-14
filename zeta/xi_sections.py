"""Hermite-function sections of the Riemann Xi function (Turan's setting).

F(t) = Xi(SCALE * t),  Xi(T) = xi(1/2 + iT)  (real, even, zeros = Riemann zeros).
Expand F = sum_k b_k psi_k in orthonormal Hermite functions (even k only).
Sections S_N = sum_{k<=N} b_k psi_k = (polynomial) * e^{-t^2/2}.
Study: coefficient decay, capture of Riemann zeros, spurious zero geometry.
"""
import mpmath as mp
import pickle, os, sys

mp.mp.dps = 80
SCALE = 10
NMAX = 60
CACHE = "xi_coeffs.pkl"


def Xi(T):
    s = mp.mpf(1)/2 + 1j*mp.mpf(T)
    return (s*(s-1)/2 * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)).real


def F(t):
    return Xi(SCALE * t)


def psi(k, t):
    # orthonormal Hermite function, stable recurrence on normalized functions:
    # psi_0 = pi^{-1/4} e^{-t^2/2};  psi_{k+1} = t*sqrt(2/(k+1)) psi_k - sqrt(k/(k+1)) psi_{k-1}
    p0 = mp.pi**mp.mpf(-0.25) * mp.e**(-t*t/2)
    if k == 0:
        return p0
    p1 = mp.sqrt(2) * t * p0
    for j in range(1, k):
        p0, p1 = p1, t*mp.sqrt(mp.mpf(2)/(j+1))*p1 - mp.sqrt(mp.mpf(j)/(j+1))*p0
    return p1


def coefficients():
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            return pickle.load(f)
    bs = {}
    for k in range(0, NMAX + 1, 2):
        b = 2 * mp.quad(lambda t: F(t) * psi(k, t), [0, 2, 4, 6, 8, 11])
        bs[k] = b
        print("b_%d = %s" % (k, mp.nstr(b, 8)), flush=True)
    with open(CACHE, "wb") as f:
        pickle.dump(bs, f)
    return bs


def hermite_coeffs(n):
    a, b = [1], [0, 2]
    if n == 0:
        return a
    for i in range(1, n):
        nxt = [0]*(i+2)
        for j, cj in enumerate(b):
            nxt[j+1] += 2*cj
        for j, cj in enumerate(a):
            nxt[j] -= 2*i*cj
        a, b = b, nxt
    return b


def section_poly(bs, N):
    """Monomial coefficients of Q_N where S_N = Q_N(t) e^{-t^2/2}."""
    coeffs = [mp.mpf(0)] * (N + 1)
    for k in range(0, N + 1, 2):
        ck = (mp.mpf(2)**k * mp.factorial(k) * mp.sqrt(mp.pi))**mp.mpf(-0.5)
        w = bs[k] * ck
        for j, hj in enumerate(hermite_coeffs(k)):
            coeffs[j] += w * hj
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


if __name__ == "__main__":
    bs = coefficients()

    # 1. decay: ln|b_k| against sqrt(2k) -- Hille-type prediction slope -pi*SCALE/4
    print("\ndecay check (prediction: slope = -pi*SCALE/4 = %.4f):" % (-mp.pi*SCALE/4))
    pts = [(mp.sqrt(2*k), mp.log(abs(bs[k]))) for k in range(8, NMAX+1, 2) if bs[k] != 0]
    for i in range(1, len(pts)):
        x0, y0 = pts[i-1]; x1, y1 = pts[i]
        if i % 5 == 0:
            print("  k=%2d  sqrt(2k)=%6.3f  ln|b_k|=%9.3f  local slope=%7.3f"
                  % ((i+4)*2, float(x1), float(y1), float((y1-y0)/(x1-x0))))

    # 2. zeros of sections
    mp.mp.dps = 120
    for N in [20, 40, 60]:
        coeffs = section_poly(bs, N)
        roots = mp.polyroots([c for c in coeffs[::-1]], maxsteps=300, extraprec=400)
        rr = sorted([float(z.real) for z in roots if abs(z.imag) < 1e-6 and z.real > 0])
        cc = [(float(z.real), float(z.imag)) for z in roots if z.imag > 1e-6]
        print("\nN=%d:  positive real zeros * SCALE:" % N)
        print("  ", ", ".join("%.4f" % (r*SCALE) for r in rr))
        print("   complex zeros (upper half, scaled by 1/sqrt(2N)):")
        sc = 1/(2*N)**0.5
        print("  ", ", ".join("(%.3f,%.3f)" % (x*sc, y*sc) for x, y in sorted(cc)))

    # true Riemann zeros for comparison
    print("\nfirst Riemann zeros:", ", ".join(mp.nstr(mp.zetazero(n).imag, 6) for n in range(1, 13)))
