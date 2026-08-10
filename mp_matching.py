"""Polynomial approximation of Meixner--Pollaczek FUNCTIONS by MP polynomials:
the MP analogue of the Hermite matching matrix.

Weight   W(x) = |Gamma(1/4 + ix)|^2   on R  (MP lambda=1/4, phi=pi/2).
MP polys P_n:  (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1},  P_0 = 1, P_1 = 2x.
MP functions   Phi_m = P_m sqrt(W) / sqrt(n_m)   (orthonormal in L^2(dx)).

Matching matrix (orthonormal form):
    b_{mk} = <Phi_m, Phat_k>_W = (n_m n_k)^{-1/2} \int P_m P_k W^{3/2} dx
           = \int Phat_m Phat_k |Gamma(1/4+ix)|^3 dx
-- the Gamma-cubed analogue of the Hermite matrix's int H_m H_k e^{-3t^2/2}.

Hermite analogue results to compare against: ratio rho = 3^{-1/2}, real spurious
band onset at |t| = sqrt((N/2) ln 3) ~ sqrt(N).  Here the weight decays like
e^{-pi|x|}, so the balance law  W(x) = rho^{2N}  predicts a LINEAR onset
x* = (2N/pi) ln(1/rho).
"""
import mpmath as mp
from fractions import Fraction
import pickle, os

mp.mp.dps = 60
KMAX = 40
CACHE = "mp_matching_coeffs.pkl"


def W(x):
    g = mp.gamma(mp.mpf(1)/4 + 1j*x)
    return (g * mp.conj(g)).real


def mp_poly_coeffs(nmax):
    """Monomial coefficients (Fractions) of P_0..P_nmax."""
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        # (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1}
        pn, pm = polys[n], polys[n-1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j+1] += 2 * c
        for j, c in enumerate(pm):
            nxt[j] -= (Fraction(n) - Fraction(1, 2)) * c
        polys.append([c / (n + 1) for c in nxt])
    return polys


POLYS = mp_poly_coeffs(KMAX)


def evalP(n, x):
    p0, p1 = mp.mpf(1), 2*x
    if n == 0:
        return p0
    for j in range(1, n):
        p0, p1 = p1, (2*x*p1 - (j - mp.mpf(1)/2)*p0) / (j + 1)
    return p1


NODES = list(range(0, 31, 2)) + [35, 40, 50]


def quadW(f):
    """integral over R of even integrand f(x) (weight included by caller)."""
    return 2 * mp.quad(f, NODES)


def norms(kmax):
    return [quadW(lambda x, n=n: evalP(n, x)**2 * W(x)) for n in range(kmax+1)]


if __name__ == "__main__":
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            nrm, B = pickle.load(f)
    else:
        nrm = norms(KMAX)
        # sanity: h_n = Gamma(n+1/2)/n! * const  -> ratio check
        print("norm ratio check (predict (n-1/2)/n):",
              mp.nstr(nrm[2]/nrm[1] / (mp.mpf(3)/4), 6), "~ 1 ?")
        B = {}
        for m in [0, 4]:
            for k in range(m % 2, KMAX + 1, 2):
                v = quadW(lambda x, m=m, k=k:
                          evalP(m, x) * evalP(k, x) * W(x)**mp.mpf(1.5))
                B[(m, k)] = v / mp.sqrt(nrm[m] * nrm[k])
                print("b_%d,%d = %s" % (m, k, mp.nstr(B[(m, k)], 8)), flush=True)
        with open(CACHE, "wb") as f:
            pickle.dump((nrm, B), f)

    print("\ndecay of b_0k (Hermite analogue: slope per k = -ln(3)/2 = -0.5493):")
    ks = sorted(k for (m, k) in B if m == 0)
    for i in range(2, len(ks), 4):
        k0, k1 = ks[i-1], ks[i]
        sl = (mp.log(abs(B[(0, k1)])) - mp.log(abs(B[(0, k0)]))) / (k1 - k0)
        print("  k=%2d  ln|b_0k|=%8.3f  slope=%8.4f  rho_local=%.4f"
              % (k1, float(mp.log(abs(B[(0, k1)]))), float(sl), float(mp.e**sl)))

    # zeros of approximants  A_N = sum_{k<=N} b_mk Phat_k  (polynomial part)
    mp.mp.dps = 90
    for m in [0, 4]:
        print("\n=== target Phi_%d ===" % m)
        if m == 4:
            r4 = mp.polyroots([mp.mpf(c.numerator)/mp.mpf(c.denominator)
                               for c in POLYS[4][::-1]])
            print("genuine zeros (P_4): "
                  + ", ".join(mp.nstr(z.real, 6) for z in r4 if z.real > 0))
        for N in [16, 24, 32, 40]:
            coeffs = [mp.mpf(0)] * (N + 1)
            for k in range(m % 2, N + 1, 2):
                wk = B[(m, k)] / mp.sqrt(nrm[k])
                for j, c in enumerate(POLYS[k]):
                    coeffs[j] += wk * mp.mpf(c.numerator) / mp.mpf(c.denominator)
            while len(coeffs) > 1 and coeffs[-1] == 0:
                coeffs.pop()
            roots = mp.polyroots(coeffs[::-1], maxsteps=300, extraprec=300)
            rr = sorted(float(z.real) for z in roots
                        if abs(z.imag) < 1e-6 and z.real > 1e-9)
            cc = sorted((float(z.real), float(z.imag))
                        for z in roots if z.imag > 1e-6)
            # balance prediction: W(x*) = b_N^2  =>  ln W(x*) = 2 ln|b_N|
            lnb = mp.log(abs(B[(m, N)]))
            xstar = mp.findroot(lambda x: mp.log(W(x)) - 2*lnb, abs(lnb))
            print("N=%d: predicted onset x* = %.2f" % (N, float(xstar)))
            print("  real zeros > 0: "
                  + ", ".join("%.4f" % r for r in rr))
            print("  complex (upper): "
                  + ", ".join("(%.2f,%.2f)" % z for z in cc))
