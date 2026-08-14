"""Dual matching matrix for Meixner-Pollaczek lambda=3/4, phi=pi/2.

Weight  W(x) = Gamma(3/4 + ix) Gamma(3/4 - ix) = |Gamma(3/4+ix)|^2.
Recurrence: (n+1) P_{n+1} = 2x P_n - (n + 1/2) P_{n-1}, P_0=1, P_1=2x.

Dual pair:
    F_m(x) = Gamma(3/4+ix) Phat_m(x)
    G_m(x) = Gamma(3/4-ix) Phat_m(x)

Dual matching coefficients:
    c_mk = int Phat_m Phat_k Gamma(3/4+ix)^2 Gamma(3/4-ix) dx

Base integral (3-Gamma lemma, a=b=c=3/4):
    J_00 = int Gamma(3/4+ix)^2 Gamma(3/4-ix) dx
         = 2pi Gamma(3/2)^2 U(3/2, 1, 1)

Now U(3/2, 1, z) = e^z/sqrt(pi) * integral form... but let's just compute it.
Using U(a, 1, z) relations: U(3/2, 1, 1) can be related to Bessel K via
U(nu+1/2, 2nu+1, 2z) = K_nu(z) / (sqrt(pi) * (2z)^nu * e^{-z})
=> with nu=1, 2nu+1=3, 2nu+1/2=2.5... doesn't fit directly.
Instead use: U(3/2, 1, 1) numerically, and the derivative algebra
dU/dz = -U(a+1, b+1, z) shifts a by 1, so the two generators are
U(3/2, 1, 1) and U(5/2, 2, 1) (or equivalently two Bessel-type constants).

Actually: U(3/2, 1, z) = e^z * Gamma(1/2, z) (since U(a,1,z) = e^z Gamma(0,z)
for a=1, and U(3/2,1,z) ~ integral...). Let's just use mpmath.

Note: U(a,b,z) in mpmath is mpmath.hyperu(a,b,z).
"""
import mpmath as mp
import pickle
import os
from fractions import Fraction

mp.mp.dps = 40

MMAX = 11
CACHE = "mp34_dual_matching.pkl"

# Integration nodes: weight decays as e^{-pi|x|/2} for lambda=3/4
# (same exponential rate as lambda=1/4, both MP with phi=pi/2)
NODES = [-50, -40, -30, -20, -12, -8, -5, -3, -1.5, 0,
         1.5, 3, 5, 8, 12, 20, 30, 40, 50]


def W(x):
    """Weight |Gamma(3/4+ix)|^2."""
    g = mp.gamma(mp.mpf(3)/4 + 1j*x)
    return (g * mp.conj(g)).real


def K34(x):
    """Three-Gamma kernel Gamma(3/4+ix)^2 Gamma(3/4-ix)."""
    return mp.gamma(mp.mpf(3)/4 + 1j*x)**2 * mp.gamma(mp.mpf(3)/4 - 1j*x)


def evalP34(n, x):
    """Evaluate P_n for lambda=3/4 at x.
    Recurrence: (n+1) P_{n+1} = 2x P_n - (n + 1/2) P_{n-1}.
    """
    if n == 0:
        return mp.mpf(1)
    p0, p1 = mp.mpf(1), 2*x
    for j in range(1, n):
        # (j+1) P_{j+1} = 2x P_j - (j + 1/2) P_{j-1}
        p0, p1 = p1, (2*x*p1 - (j + mp.mpf(1)/2)*p0) / (j + 1)
    return p1


def mp34_poly_coeffs(nmax):
    """Monomial coefficients (Fractions) of P_0..P_nmax for lambda=3/4."""
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        # (n+1) P_{n+1} = 2x P_n - (n + 1/2) P_{n-1}
        pn, pm = polys[n], polys[n-1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j+1] += 2 * c
        for j, c in enumerate(pm):
            nxt[j] -= (Fraction(n) + Fraction(1, 2)) * c
        polys.append([c / (n + 1) for c in nxt])
    return polys


POLYS34 = mp34_poly_coeffs(MMAX + 5)


def main():
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            nrm, C = pickle.load(f)
        print(f"Loaded cache: {CACHE}")
    else:
        print("Computing norms h_n = int P_n^2 W dx ...")
        nrm = []
        for n in range(MMAX + 1):
            # Weight decays, use even+odd quadrature
            val = mp.quad(lambda x, n=n: evalP34(n, x)**2 * W(x), NODES)
            nrm.append(val)
            print(f"  h_{n} = {mp.nstr(val, 12)}")

        print("\nComputing matching matrix c_mk ...")
        C = {}
        for m in range(MMAX + 1):
            for k in range(m, MMAX + 1):
                J = mp.quad(lambda x, m=m, k=k:
                            evalP34(m, x) * evalP34(k, x) * K34(x), NODES)
                cmk = J / mp.sqrt(nrm[m] * nrm[k])
                C[(m, k)] = cmk
                print(f"  c_{m},{k} = {mp.nstr(cmk, 10)}", flush=True)

        with open(CACHE, "wb") as f:
            pickle.dump((nrm, C), f)
        print(f"\nSaved to {CACHE}")

    def c(m, k):
        return C[(min(m, k), max(m, k))]

    # ---- Sanity 1: orthonormality ----
    print("\n=== Orthonormality check ===")
    # Check h_n formula: for lambda, h_n = Gamma(n+2*lambda)/(n! * Gamma(2*lambda)) * h_0
    # For lambda=3/4: h_n = Gamma(n+3/2)/(n! * Gamma(3/2)) * h_0
    # = sqrt(pi)/2 * Gamma(n+3/2) / (n! * sqrt(pi)/2) = Gamma(n+3/2)/n!
    # But wait, we need to check the actual formula.
    # Orthogonality weight integral: int P_m P_k W dx = h_n delta_{mk}
    # for MP: h_n = (2pi)/(2lambda)_n / n! * Gamma(2lambda) ... let's check ratio
    for n in range(1, 5):
        ratio = nrm[n] / nrm[n-1]
        # Expected: h_n/h_{n-1} = (n - 1 + 2*lambda)/n = (n + 1/2)/n for lambda=3/4
        expected = (n + mp.mpf(1)/2) / n
        print(f"  h_{n}/h_{n-1} = {mp.nstr(ratio, 10)} (expected {mp.nstr(expected, 8)})")

    # ---- Sanity 2: base integral vs 3-Gamma lemma ----
    print("\n=== Base integral J_00 ===")
    J00_numerical = C[(0, 0)] * nrm[0]
    # 3-Gamma lemma: a=b=c=3/4
    # J_00 = 2pi * Gamma(3/2)^2 * U(3/2, 1, 1)
    J00_closed = 2 * mp.pi * mp.gamma(mp.mpf(3)/2)**2 * mp.hyperu(mp.mpf(3)/2, 1, 1)
    print(f"  Numerical:  {mp.nstr(J00_numerical, 20)}")
    print(f"  3-Gamma lemma: {mp.nstr(J00_closed, 20)}")
    print(f"  Difference: {mp.nstr(abs(J00_numerical - J00_closed), 5)}")

    # Also try expressing via Bessel K
    # U(3/2, 1, z): use mpmath to evaluate and try to express via Bessel
    U_val = mp.hyperu(mp.mpf(3)/2, 1, 1)
    print(f"\n  U(3/2, 1, 1) = {mp.nstr(U_val, 20)}")
    # U(3/2, 2, 2z) = K_1(z) / (sqrt(pi) * e^{-z} * (2z))... let's try
    # U(nu+1/2, 2nu+1, 2z) = K_nu(z)/(sqrt(pi)*(2z)^nu * e^z)
    # nu+1/2 = 3/2 => nu=1; 2nu+1=3 != 1. So U(3/2, 1, z) is not directly K_1.
    # Try U(3/2, 1, 1) numerically vs Gamma/E1 type
    # U(a, 1, z) = -1/Gamma(a) * (ln(z)*M(a,1,z) + F(a,1,z)) -- complex
    # Let's just check if it's related to K_0, K_1
    A_prime = mp.sqrt(2*mp.e) * mp.besselk(0, mp.mpf(1)/2)
    B_prime = mp.sqrt(2*mp.e) * mp.besselk(1, mp.mpf(1)/2)
    print(f"\n  A' = sqrt(2e)*K_0(1/2) = {mp.nstr(A_prime, 15)}")
    print(f"  B' = sqrt(2e)*K_1(1/2) = {mp.nstr(B_prime, 15)}")
    # Try: U(3/2,1,1) = a*K_0(1/2) + b*K_1(1/2) via PSLQ
    pslq_U = mp.pslq([U_val, A_prime, B_prime], tol=1e-15, maxcoeff=100)
    print(f"  PSLQ U(3/2,1,1) vs K_0(1/2),K_1(1/2): {pslq_U}")

    # ---- checkerboard / phase law ----
    print("\n=== Checkerboard check ===")
    worst_cb = max(abs(mp.im(c(m, k))) if (m + k) % 2 == 0 else abs(mp.re(c(m, k)))
                   for m in range(MMAX + 1) for k in range(MMAX + 1))
    print(f"  Max violation (should be ~0): {mp.nstr(worst_cb, 4)}")

    print("\n=== Phase law check: c_mk = (-i)^{|m-k|} |c_mk| ? ===")
    worst_phase = mp.mpf(0)
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            diff = abs(m - k)
            expected_phase = (-1j)**diff
            v = c(m, k)
            if abs(v) > 1e-30:
                # expected_phase * |c_mk| should = c_mk
                reconstructed = expected_phase * abs(v)
                err = abs(v - reconstructed)
                worst_phase = max(worst_phase, err)
    print(f"  Max phase-law violation: {mp.nstr(worst_phase, 4)}")

    # ---- the table: raw c_mk ----
    print("\n=== c_mk table (m+k even: real; m+k odd: imaginary) ===")
    print("     " + "".join(f"k={k:<13d}" for k in range(min(7, MMAX + 1))))
    for m in range(min(7, MMAX + 1)):
        row = []
        for k in range(min(7, MMAX + 1)):
            v = c(m, k)
            if (m + k) % 2 == 0:
                row.append(f"{float(mp.re(v)):+.5f}      ")
            else:
                row.append(f"{float(mp.im(v)):+.5f}i     ")
        print(f"m={m:<3d} " + "".join(row))

    # ---- decay of row 0 ----
    print("\n=== Decay of |c_0k| ===")
    for k in range(0, MMAX + 1):
        v = abs(c(0, k))
        print(f"  k={k:2d}  |c_0k| = {mp.nstr(v, 10)}")

    # ---- PSLQ for integrality ----
    print("\n=== PSLQ integrality check ===")
    h0 = nrm[0]
    # For lambda=1/4: T_mk = 2^m m! 2^k k! / h0 * J_mk
    # The J_mk = c_mk * sqrt(h_m h_k)
    # For lambda=3/4, try same normalization first
    # h_n for lambda=3/4: h_n = h_0 * prod_{j=1}^n (j+1/2)/j = h_0 * Gamma(n+3/2)/(n! * Gamma(3/2))
    # = h_0 * (2n+1)!! / (2^n * n! * sqrt(pi)/2 * 2^n/sqrt(pi)) ... let's compute
    def h_formula(n):
        # h_n = h_0 * C(2n,n) * (2/4)^n ... for lambda=1/4 it was sqrt(2)*pi^{3/2} * C(2n,n)/4^n
        # For lambda=3/4: different formula
        # h_n = (2*lambda)_n / n! * h_0 * something... use the ratio
        return nrm[n]

    # Candidates for integrality basis
    # For lambda=3/4: U(3/2,1,1) and its contiguous U(5/2,2,1)?
    # Or expressed via Bessel: U(a,1,z) relates to K_0, E_1...
    # Let's try K_{1/4} type... or just try numeric PSLQ
    U32 = mp.hyperu(mp.mpf(3)/2, 1, 1)
    U52 = mp.hyperu(mp.mpf(5)/2, 2, 1)
    # Also try U(3/2,1,1) expressed via K_1 with half-integer args
    # K_{1/2}(z) = sqrt(pi/(2z)) e^{-z}, elementary
    # K_{3/2}(z) = sqrt(pi/(2z)) e^{-z}(1 + 1/z)
    Khalf = mp.besselk(mp.mpf(1)/2, mp.mpf(1)/2)
    K3half = mp.besselk(mp.mpf(3)/2, mp.mpf(1)/2)
    print(f"  U(3/2, 1, 1) = {mp.nstr(U32, 15)}")
    print(f"  U(5/2, 2, 1) = {mp.nstr(U52, 15)}")
    print(f"  K_{{1/2}}(1/2) = {mp.nstr(Khalf, 15)}")
    print(f"  K_{{3/2}}(1/2) = {mp.nstr(K3half, 15)}")

    # For each (m,k), compute T_mk = 2^m m! 2^k k! / h0 * J_mk
    # where J_mk = c_mk * sqrt(h_m * h_k)
    print("\n  Testing T_mk = 2^m m! 2^k k! / h0 * J_mk ∈ Z*A' + Z*B'")
    print("  (phase removed: T_mk / (-i)^{|m-k|} should be real)")

    # First identify A', B' for lambda=3/4
    # Try: J_00 = 2pi Gamma(3/2)^2 U(3/2,1,1)
    # = 2pi * (sqrt(pi)/2)^2 * U(3/2,1,1) = pi^2/2 * U(3/2,1,1)
    # So h0 = J_00 / c00 ... wait c00 is normalized so J_00 = c00 * h0
    # T_00 = 1! * 1! / h0 * J_00 = J_00/h0 = c00
    # Hmm, let's think more carefully.
    # For lambda=1/4: T_mk = 2^m m! 2^k k! / h0 * J_mk, and T_00 = A = sqrt(2e)*K_0(1/2)
    # Because J_00 / h0 = c00 = 1 (by construction of Phat? No: c00 = J_00/(h0) since
    # c_mk = J_mk / sqrt(h_m * h_k), so c_00 = J_00/h_0... yes because sqrt(h_0*h_0) = h_0.
    # So T_00 = 1*1/h_0 * J_00 = c_00. Let me verify by checking the lambda=1/4 code.
    # In mp_dual_matching.py: C[(m,k)] = J / sqrt(nrm[m]*nrm[k]). So c_00 = J_00/h_0.
    # T_00 = 2^0*0! * 2^0*0! / h0 * J_00 = J_00/h_0 = c_00.
    # For lambda=3/4 likewise, T_00 = c_00.

    # Let's try different normalizations to find integrality
    # Conjecture: same structure but different A', B'
    # A' = U(3/2,1,1) times some normalization
    # From J_00 = 2pi*Gamma(3/2)^2*U(3/2,1,1):
    # c_00 = J_00/h_0.
    # T_00 = c_00 = J_00/h_0.
    # For the denominator: the recurrence for lambda=3/4 has (n+1/2) instead of (n-1/2).
    # h_n/h_0 for lambda=3/4 = prod_{j=1}^n (j+1/2)/j = (3/2)(5/2)...(n+1/2) / n!
    # = Gamma(n+3/2) / (Gamma(3/2) * n!)
    # = (2n+1)!! / (2^n * n! * sqrt(pi)/2) ... let me just compute
    h_ratios = [nrm[n] / nrm[0] for n in range(MMAX + 1)]
    print("\n  h_n / h_0 ratios:")
    for n in range(6):
        print(f"    h_{n}/h_0 = {mp.nstr(h_ratios[n], 10)}")
        # Compare with Gamma(n+3/2)/(n! * Gamma(3/2))
        expected_ratio = mp.gamma(n + mp.mpf(3)/2) / (mp.fac(n) * mp.gamma(mp.mpf(3)/2))
        print(f"    Gamma(n+3/2)/(n! Gamma(3/2)) = {mp.nstr(expected_ratio, 10)}")

    # So h_n = h_0 * Gamma(n+3/2) / (n! * Gamma(3/2)) = h_0 * (2n+1)!! / (2^n * n! * (1/2))
    # T_mk with same denominator as lambda=1/4: 2^m m! * 2^k k!
    # But h_n for lambda=3/4 scales differently. Let's try:
    # T_mk^(3/4) = 2^m m! * 2^k k! / h_0 * J_mk
    # = 2^m m! * 2^k k! * c_mk * sqrt(h_m h_k) / h_0
    # = 2^m m! * 2^k k! * c_mk * sqrt(h_m/h_0 * h_k/h_0)
    print("\n  Computing T_mk = 2^m m! 2^k k! * c_mk * sqrt(h_m*h_k) / h_0")
    print("  and running PSLQ against U(3/2,1,1), U(5/2,2,1):")

    # Define candidate constants
    # The derivative of U(a,b,z) = -a*U(a+1,b+1,z)
    # so dU/dz(3/2,1,z)|_{z=1} = -3/2 * U(5/2,2,1)
    # The GF will be e^{zeta/2} U(3/2,1,zeta) type -- need two generators
    # Let's pick A' = 2pi*Gamma(3/2)^2 * U(3/2,1,1) / h_0 = J_00/h_0 = c_00
    # and B' from the contiguous: d/dz U(3/2,1,z) = -3/2 * U(5/2,2,1)
    # But we need B' to be a "nice" constant. Let's see what PSLQ finds.

    AP = C[(0, 0)]   # = J_00/h_0, our A' candidate (like A for lambda=1/4)
    print(f"\n  A' candidate = c_00 = {mp.nstr(AP, 15)}")

    # For each entry, compute T_mk (phase-corrected, real)
    pslq_hits = {}
    for m in range(6):
        for k in range(m, 6):
            J_mk = C[(m, k)] * mp.sqrt(nrm[m] * nrm[k])
            T_raw = (2**m * mp.fac(m) * 2**k * mp.fac(k)) / nrm[0] * J_mk
            # Remove phase
            diff = abs(m - k)
            phase = (-1j)**diff
            T_real = T_raw / phase
            # Should be real
            if abs(mp.im(T_real)) > 1e-20:
                print(f"  WARNING: T_{m},{k} not real! Im = {mp.nstr(mp.im(T_real), 5)}")
            T = mp.re(T_real)
            # Try PSLQ vs A' and various U values
            U32_val = mp.hyperu(mp.mpf(3)/2, 1, 1)
            U52_val = mp.hyperu(mp.mpf(5)/2, 2, 1)
            rel = mp.pslq([T, U32_val, U52_val], tol=1e-20, maxcoeff=10000)
            pslq_hits[(m, k)] = (T, rel)
            rel_str = str(rel) if rel else "not found"
            print(f"  T_{m},{k} = {mp.nstr(T, 10)}  PSLQ: {rel_str}")

    return nrm, C, pslq_hits


if __name__ == "__main__":
    nrm, C, pslq_hits = main()
