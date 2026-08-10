"""Polynomial approximation (log-Hermite polynomials) of the Riemann Xi function.

Contrast with xi_sections.py (basis-expansion sections S_N = poly * e^{-t^2/2}):
here we take the best L^2(R, e^{-t^2} dt)-approximation of F(t) = Xi(SCALE*t)
by polynomials of degree <= N.  On R_+ under x = e^t this is approximation by
polynomials in log x against the log-normal weight ("log-Hermite polynomials").

Key identity: with orthonormal Hermite polys  h_k = H_k / (2^k k! sqrt(pi))^{1/2}
= psi_k * e^{t^2/2},  the projection coefficients are

    beta_k = <F, h_k>_{e^{-t^2}} = \int F(t) psi_k(t) e^{-t^2/2} dt,

i.e. the Hermite-FUNCTION coefficients of the damped target F * e^{-t^2/2}.
Prediction (zero.tex taxonomy): beta_k ~ 3^{-k/2} (matching-matrix universality),
real spurious zeros past the balance height  t^2/2 + (pi*SCALE/4) t = (N/2) ln 3.
"""
import mpmath as mp
import pickle, os

mp.mp.dps = 80
SCALE = 10
NMAX = 60
CACHE = "xi_polyapprox_coeffs.pkl"


def Xi(T):
    s = mp.mpf(1)/2 + 1j*mp.mpf(T)
    return (s*(s-1)/2 * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)).real


def F(t):
    return Xi(SCALE * t)


def psi(k, t):
    # orthonormal Hermite function, stable normalized recurrence
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
        # beta_k = int F(t) psi_k(t) e^{-t^2/2} dt   (even integrand)
        b = 2 * mp.quad(lambda t: F(t) * psi(k, t) * mp.e**(-t*t/2),
                        [0, 2, 4, 6, 8, 11])
        bs[k] = b
        print("beta_%d = %s" % (k, mp.nstr(b, 8)), flush=True)
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


def approx_poly(bs, N):
    """Monomial coefficients of P_N F = sum_{k<=N} beta_k h_k (a pure polynomial)."""
    coeffs = [mp.mpf(0)] * (N + 1)
    for k in range(0, N + 1, 2):
        ck = (mp.mpf(2)**k * mp.factorial(k) * mp.sqrt(mp.pi))**mp.mpf(-0.5)
        w = bs[k] * ck
        for j, hj in enumerate(hermite_coeffs(k)):
            coeffs[j] += w * hj
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def predicted_onset(bs, N):
    """Solve t^2/2 + (pi*SCALE/4) t = -ln|beta_N| (crude balance, no poly factors)."""
    rhs = -mp.log(abs(bs[N]))
    a, b = mp.mpf(1)/2, mp.pi*SCALE/4
    t = (-b + mp.sqrt(b*b + 4*a*rhs)) / (2*a)
    return t


if __name__ == "__main__":
    bs = coefficients()

    # 1. decay: prediction ln|beta_k| ~ -(k/2) ln 3, slope per k = -0.5493
    print("\ndecay check (prediction: slope per k -> -ln(3)/2 = %.4f):"
          % float(-mp.log(3)/2))
    ks = [k for k in range(0, NMAX+1, 2) if bs[k] != 0]
    for i in range(2, len(ks), 5):
        k0, k1 = ks[i-1], ks[i]
        sl = (mp.log(abs(bs[k1])) - mp.log(abs(bs[k0]))) / (k1 - k0)
        print("  k=%2d  ln|beta_k|=%9.3f  local slope=%7.4f"
              % (k1, float(mp.log(abs(bs[k1]))), float(sl)))

    # 2. zeros of polynomial approximants
    mp.mp.dps = 120
    for N in [20, 40, 60]:
        coeffs = approx_poly(bs, N)
        roots = mp.polyroots([c for c in coeffs[::-1]], maxsteps=300, extraprec=400)
        rr = sorted([float(z.real) for z in roots if abs(z.imag) < 1e-6 and z.real > 0])
        cc = [(float(z.real), float(z.imag)) for z in roots if z.imag > 1e-6]
        tc = predicted_onset(bs, N)
        print("\nN=%d:  predicted onset T_c ~ %.1f  (t_c=%.2f, turning point T=%.0f)"
              % (N, float(tc*SCALE), float(tc), float(SCALE*mp.sqrt(2*N))))
        print("   positive real zeros * SCALE:")
        print("  ", ", ".join("%.4f" % (r*SCALE) for r in rr))
        print("   complex zeros (upper half, * SCALE):")
        print("  ", ", ".join("(%.2f,%.2f)" % (x*SCALE, y*SCALE) for x, y in sorted(cc)))

    print("\nfirst Riemann zeros:",
          ", ".join(mp.nstr(mp.zetazero(n).imag, 6) for n in range(1, 13)))
