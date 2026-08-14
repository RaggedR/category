"""
Jacobi matching matrix — endpoint-singularity control case.

Weight W_αβ(x) = (1−x)^α (1+x)^β on [−1,1].
Orthonormal Jacobi polynomials P̂ₙ^{α,β}.
Matching matrix b_mk = ∫_{-1}^{1} P̂ₘ P̂ₖ W^{3/2} dx.

Four cases studied:
  (A) Chebyshev-T: α=β=−½  (W^{3/2} = (1−x²)^{−3/4})
  (B) Chebyshev-U: α=β=+½  (W^{3/2} = (1−x²)^{+3/4})
  (C) Legendre:    α=β=0   (degenerate: W^{3/2}=W=1, B=I)
  (D) Asymmetric:  α=½, β=−½

Uses x=cosθ substitution for Chebyshev cases to get Γ-function closed forms.
"""

import mpmath as mp
import numpy as np
import pickle, os, sys

mp.mp.dps = 50
NMAX = 11   # 12×12 table (indices 0..11)
CACHE = "/Users/robin/git/category/jacobi_matching.pkl"

# ─── Jacobi polynomial recurrence ──────────────────────────────────────────
def jacobi_poly(n, a, b, x):
    """Evaluate P_n^{(a,b)}(x) via 3-term recurrence."""
    if n == 0:
        return mp.mpf(1)
    p0 = mp.mpf(1)
    p1 = (a - b)/2 + (a + b + 2)/2 * x
    for k in range(1, n):
        A = 2*(k+1)*(k+a+b+1)*(2*k+a+b)
        B = (2*k+a+b+1)*(a**2 - b**2)
        C = (2*k+a+b)*(2*k+a+b+1)*(2*k+a+b+2)
        D = 2*(k+a)*(k+b)*(2*k+a+b+2)
        p0, p1 = p1, ((B + C*x)*p1 - D*p0) / A
    return p1

def jacobi_norm_sq(n, a, b):
    """h_n = ||P_n^{(a,b)}||^2 = ∫_{-1}^1 [P_n]^2 (1-x)^a (1+x)^b dx."""
    from mpmath import gamma, factorial, power
    if n == 0:
        return mp.power(2, a+b+1) * mp.gamma(a+1)*mp.gamma(b+1)/mp.gamma(a+b+2)
    num = mp.power(2, a+b+1) * mp.gamma(n+a+1)*mp.gamma(n+b+1)
    den = mp.factorial(n) * (2*n+a+b+1) * mp.gamma(n+a+b+1)
    return num / den

def jacobi_normalized(n, a, b, x):
    """P̂_n(x) = P_n(x) / sqrt(h_n)."""
    h = jacobi_norm_sq(n, a, b)
    return jacobi_poly(n, a, b, x) / mp.sqrt(h)

# ─── Numerical integration of b_mk ─────────────────────────────────────────
def compute_b_numerical(m, k, alpha, beta):
    """b_mk = ∫_{-1}^1 P̂ₘ^{α,β}(x) P̂ₖ^{α,β}(x) (1-x)^{3α/2} (1+x)^{3β/2} dx"""
    a3, b3 = 3*alpha/2, 3*beta/2
    def integrand(x):
        Pm = jacobi_normalized(m, alpha, beta, x)
        Pk = jacobi_normalized(k, alpha, beta, x)
        # weight W^{3/2} = (1-x)^{3α/2}(1+x)^{3β/2}
        if x >= mp.mpf('1') or x <= mp.mpf('-1'):
            return mp.mpf(0)
        wt = mp.power(1 - x, a3) * mp.power(1 + x, b3)
        return Pm * Pk * wt
    # Gauss-Legendre / tanh-sinh works well for integrable singularities
    val = mp.quad(integrand, [-1, 1], error=False, maxdegree=8)
    return val

# ─── CHEBYSHEV-T CLOSED FORMS ───────────────────────────────────────────────
# W = (1-x²)^{-1/2}, normalized polys: P̂₀ = 1/√π, P̂ₙ = √(2/π) Tₙ(x) n≥1
# x = cosθ: ∫_{-1}^1 T̂ₘ T̂ₖ (1-x²)^{-3/4} dx
#           = (norm factors) ∫₀^π cos(mθ) cos(kθ) sin^{-1/2}(θ) dθ
# The integral I(m,k,p) = ∫₀^π cos(mθ)cos(kθ) sinᵖ(θ) dθ
# = (π/2) · 2^{-p} · B((p+1)/2, 1/2) · [... trigonometric convolution formula]
# Exact: I(m,k,p) = (π/2^{p+1}) · Γ(p+1)/[Γ((p+m+k+2)/2)Γ((p-m-k+2)/2)]
#   ... uses the formula for ∫₀^π e^{i(m-k)θ} sinᵖ(θ) dθ.
#
# More precisely, ∫₀^π cos(mθ)cos(kθ) sinᵖ(θ) dθ
# = (π/4) · [B((p+1)/2, |(m-k)|+... ) + B(...)]-type expression.
# We use the Wallis-type formula directly:
#   ∫₀^π sinᵖ(θ) cos(qθ) dθ = cos(πq/2)·π·Γ(p+1) / (2^p · Γ((p+q+2)/2) · Γ((p-q+2)/2))
# This is valid for p > -1, q integer ≥ 0.  (Phase factor added 2026-08-10.)
# cos(mθ)cos(kθ) = [cos((m+k)θ) + cos((m-k)θ)] / 2

def I_sin_cos(p, q):
    """∫₀^π sin^p(θ) cos(qθ) dθ  (q ≥ 0 integer, p > -1)
    = cos(πq/2) · π Γ(p+1) / (2^p Γ((p+q+2)/2) Γ((p-q+2)/2))
    CORRECTED 2026-08-10: the phase factor cos(πq/2) is essential.
    It kills all odd q and contributes (−1)^{q/2} for even q, which
    cancels the alternating sign of Γ((p−q+2)/2).
    """
    p = mp.mpf(p)
    q = int(abs(q))
    if q % 2 == 1:
        return mp.mpf(0)
    phase = (-1) ** (q // 2)
    num = phase * mp.pi * mp.gamma(p + 1)
    den = mp.power(2, p) * mp.gamma((p + q + 2)/2) * mp.gamma((p - q + 2)/2)
    return num / den

def b_chebyshev_T_closed(m, k):
    """
    Chebyshev-T: α=β=−½.
    W = (1-x²)^{-1/2}/π (the standard Chebyshev weight, but for orthonormality
    we use the usual normalizations).
    P̂₀ = 1/√π, P̂ₙ = √(2/π) Tₙ(x) for n≥1.
    b_mk = ∫_{-1}^1 P̂ₘ P̂ₖ (1-x²)^{-3/4} dx
    Under x=cosθ, dx = sinθ dθ, (1-x²)^{-3/4} = sin^{-3/2}(θ):
    b_mk = norm_m · norm_k · ∫₀^π cos(mθ)cos(kθ) sin^{-3/2}(θ)·sinθ dθ
         = norm_m · norm_k · ∫₀^π cos(mθ)cos(kθ) sin^{-1/2}(θ) dθ
    where norm_0 = 1/√π, norm_n = √(2/π) for n≥1.
    sin^{-1/2}(θ): p = -1/2 in I_sin_cos.
    ∫₀^π cos(mθ)cos(kθ) sin^{-1/2}(θ) dθ
      = [I_sin_cos(-1/2, m+k) + I_sin_cos(-1/2, |m-k|)] / 2
    """
    p = mp.mpf(-1)/2
    integral = (I_sin_cos(p, m+k) + I_sin_cos(p, abs(m-k))) / 2
    nm = 1/mp.sqrt(mp.pi) if m == 0 else mp.sqrt(2/mp.pi)
    nk = 1/mp.sqrt(mp.pi) if k == 0 else mp.sqrt(2/mp.pi)
    return nm * nk * integral

def b_chebyshev_U_closed(m, k):
    """
    Chebyshev-U: α=β=+½.
    P̂ₙ^{1/2,1/2}(x) relates to Uₙ(x) = sin((n+1)θ)/sinθ.
    Norm: h_n = π/4 for all n (from standard formula: 2^{-2}·π·1·1/((2n+3)/2·1)
    -> actually h_n = π/4 for all n when α=β=1/2... let's check:
    h_n = 2^{α+β+1} Γ(n+α+1)Γ(n+β+1) / (n! (2n+α+β+1) Γ(n+α+β+1))
    α=β=1/2: = 2^2 Γ(n+3/2)^2 / (n! (2n+2) Γ(n+2))
    Γ(n+3/2) = (n+1/2)! = Γ(3/2)·(3/2)···(n+1/2) = √π/2^n · (2n+1)!!/(2n·...)
    Actually: h_n = π/(n+1)... let me compute directly from norm_sq.
    x = cosθ: Uₙ(cosθ) = sin((n+1)θ)/sinθ.
    P_n^{1/2,1/2}(x) = (n+1)/(n+1) stuff... let's use norm_sq directly.
    P̂ₙ = Uₙ(x) / sqrt(h_n),  h_n^U = π/2.
    x=cosθ: Uₙ(cosθ)sinθ = sin((n+1)θ).
    b_mk = ∫_{-1}^1 P̂ₘ P̂ₖ (1-x²)^{3/4} dx
    x=cosθ: (1-x²)^{3/4} = sin^{3/2}(θ), dx=sinθ dθ
    = (1/h_U) ∫₀^π Uₘ(cosθ) Uₖ(cosθ) sin^{3/2+1}(θ) dθ
    Uₙ(cosθ) = sin((n+1)θ)/sinθ, so Uₘ Uₖ sin² = sin((m+1)θ)sin((k+1)θ)
    = [cos((m-k)θ) - cos((m+k+2)θ)] / 2
    So b_mk = (1/(π/2)) · (1/2) ∫₀^π [cos((m-k)θ) - cos((m+k+2)θ)] sin^{3/2}(θ) dθ
    """
    # h_n^{(1/2,1/2)} from formula
    # For Chebyshev-U: the standard orthonormal polys use h_n = π/2
    hU = mp.pi / 2
    p = mp.mpf(3)/2
    q1 = abs(m - k)
    q2 = m + k + 2
    integral = (I_sin_cos(p, q1) - I_sin_cos(p, q2)) / 2
    return integral / hU

# ─── BUILD TABLES ───────────────────────────────────────────────────────────
def build_table(alpha, beta, label, use_closed=None, nmax=NMAX):
    """Build (nmax+1)×(nmax+1) matching matrix by numerical quadrature."""
    n = nmax + 1
    B = mp.matrix(n, n)
    for m in range(n):
        for k in range(m, n):
            if use_closed is not None:
                val = use_closed(m, k)
            else:
                val = compute_b_numerical(m, k, alpha, beta)
            B[m, k] = val
            B[k, m] = val
        print(f"  {label}: row {m} done", flush=True)
    return B

def build_all_tables(nmax=NMAX):
    if os.path.exists(CACHE):
        print("Loading from cache...")
        with open(CACHE, "rb") as f:
            data = pickle.load(f)
        def list_to_mat(L):
            n = len(L)
            M = mp.matrix(n, n)
            for i in range(n):
                for j in range(n):
                    M[i,j] = mp.mpf(L[i][j])
            return M
        return {k: list_to_mat(v) for k, v in data.items()}

    print("Computing Chebyshev-T (closed form)...")
    BT = build_table(-mp.mpf(1)/2, -mp.mpf(1)/2, "ChebT",
                     use_closed=b_chebyshev_T_closed, nmax=nmax)

    print("Computing Chebyshev-U (closed form)...")
    BU = build_table(mp.mpf(1)/2, mp.mpf(1)/2, "ChebU",
                     use_closed=b_chebyshev_U_closed, nmax=nmax)

    print("Computing Legendre (α=β=0) — degenerate (B=I)...")
    n = nmax + 1
    BL = mp.matrix(n, n)
    for i in range(n):
        BL[i, i] = mp.mpf(1)
    # Note: W=1, W^{3/2}=1, b_mk=δ_mk. Document as fixed point.

    print("Computing asymmetric Jacobi α=1/2, β=-1/2 (numerical)...")
    BA = build_table(mp.mpf(1)/2, -mp.mpf(1)/2, "Asym", nmax=nmax)

    # Convert mpmath matrices to nested lists for pickling
    def mat_to_list(M):
        return [[float(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]
    result_lists = {"T": mat_to_list(BT), "U": mat_to_list(BU),
                    "L": mat_to_list(BL), "A": mat_to_list(BA)}
    with open(CACHE, "wb") as f:
        pickle.dump(result_lists, f)
    print("Saved to cache.")
    return {"T": BT, "U": BU, "L": BL, "A": BA}

# ─── DECAY ANALYSIS ─────────────────────────────────────────────────────────
def analyze_decay(B, label):
    """Extract b_{0k} and fit power law k^{-p}."""
    n = B.rows
    data = []
    for k in range(1, n):
        val = abs(B[0, k])
        if val > 1e-100:
            data.append((k, val))
    if len(data) < 3:
        print(f"  {label}: too few nonzero entries")
        return None
    # fit log|b_{0k}| = -p log(k) + const via last few points
    ks = np.array([d[0] for d in data])
    vs = np.array([float(mp.log(d[1])) for d in data])
    # linear fit in log-log
    coeffs = np.polyfit(np.log(ks[2:]), vs[2:], 1)
    p_fit = -coeffs[0]
    print(f"  {label}: power-law exponent p ≈ {p_fit:.4f}")
    return p_fit, data

# ─── GENERATOR RING ANALYSIS ────────────────────────────────────────────────
def find_generator(BT, BU):
    """
    For Chebyshev-T: b_mk = (norms) · [I(-1/2, m+k) + I(-1/2, |m-k|)] / 2
    The integral I(-1/2, q) = π Γ(1/2) / (2^{-1/2} Γ((1+q)/2) Γ((2-q)/2+...))
    = π^{3/2} / (2^{-1/2} · Γ((q+1)/2) Γ((1-q)/2))  — via Γ(1/2)=√π.
    The overall constant involves Γ(1/4) (from q=0 term: Γ(1/4) appears in
    ∫₀^π sin^{-1/2}(θ) dθ = √π Γ(1/4)/Γ(3/4) type).
    Check: b_{00} / b_{01} should be rational if single generator.
    """
    b00 = BT[0,0]
    b01 = BT[0,1]
    b11 = BT[1,1]
    print(f"  b_00 = {mp.nstr(b00, 10)}")
    print(f"  b_01 = {mp.nstr(b01, 10)}")
    print(f"  b_11 = {mp.nstr(b11, 10)}")
    # ratio b_00 / b_01
    rat = b00 / b01
    print(f"  b_00/b_01 = {mp.nstr(rat, 10)}")
    # Is b_00 a rational multiple of Γ(1/4)^2/π or similar?
    G14 = mp.gamma(mp.mpf(1)/4)
    G34 = mp.gamma(mp.mpf(3)/4)
    # γ = √π · Γ(1/4) / (2^{1/4} · Γ(3/4)·π) or similar
    # b_00 for Cheb-T: norm_0^2 * I(-1/2, 0) = (1/π) * π*Γ(1/2)/(2^{-1/2}*Γ(3/4)^2)
    # = Γ(1/2)/(2^{-1/2}*π * Γ(3/4)^2) = √π * √2 / (π * Γ(3/4)^2)
    # = √2/(√π * Γ(3/4)^2)
    b00_theory = mp.sqrt(2) / (mp.sqrt(mp.pi) * G34**2)
    print(f"  b_00 theory = {mp.nstr(b00_theory, 10)}")
    print(f"  match: {mp.nstr(b00 - b00_theory, 5)}")
    # Single generator: γ = Γ(1/4)/Γ(3/4) = π/Γ(3/4)^2 (by reflection)
    # b_{mk} = rational · γ^?
    # Let's check b_{mk} / b_{00}
    print(f"\n  b_mk / b_00 for small m,k (should be rational if one generator):")
    for m in range(4):
        for k in range(m, 4):
            ratio = BT[m,k] / b00
            print(f"    b_{m}{k}/b_00 = {mp.nstr(ratio, 8)}")

# ─── ZEROS OF TRUNCATED APPROXIMANTS ────────────────────────────────────────
def compute_zeros(BT, alpha, beta, label, nmax=NMAX):
    """
    Compute zeros of truncated expansions A_N(x) = ∑_{k=0}^N b_{0k} P̂_k(x)
    for N = 10, 20, 40 (but we only have NMAX=11, so N=10).
    """
    import numpy as np
    # We need Chebyshev expansion: A_N(x) = ∑ b_{0k} P̂_k(x)
    # For Chebyshev-T: P̂_k = (1/√π) T_k for k=0, √(2/π) T_k for k≥1
    # So A_N(x) = b_00/√π + ∑_{k=1}^N b_{0k}√(2/π) T_k(x)
    # In Chebyshev coefficients: c_0 = b_00/√π, c_k = b_{0k}√(2/π) for k≥1
    N = BT.rows - 1
    c = []
    for k in range(N+1):
        if k == 0:
            c.append(float(BT[0,0]) / np.sqrt(np.pi))
        else:
            c.append(float(BT[0,k]) * np.sqrt(2/np.pi))
    # Convert Chebyshev coefficients to monomial, then find roots
    # Or use numpy Chebyshev
    poly = np.polynomial.chebyshev.Chebyshev(c)
    # For complex zeros, we need to evaluate off [−1,1]
    companion_roots = poly.roots()
    return companion_roots

# ─── MAIN ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("JACOBI MATCHING MATRIX COMPUTATION")
    print("=" * 60)

    tables = build_all_tables(nmax=NMAX)
    BT, BU, BL, BA = tables["T"], tables["U"], tables["L"], tables["A"]

    print("\n--- Chebyshev-T (α=β=−½) first row ---")
    for k in range(NMAX+1):
        print(f"  b_0{k:2d} = {mp.nstr(BT[0,k], 12)}")

    print("\n--- Chebyshev-U (α=β=+½) first row ---")
    for k in range(NMAX+1):
        print(f"  b_0{k:2d} = {mp.nstr(BU[0,k], 12)}")

    print("\n--- Legendre (α=β=0): B = I (fixed point) ---")
    print("  b_00 = 1, b_01 = 0, ... (identity matrix)")

    print("\n--- Asymmetric (α=½, β=−½) first row ---")
    for k in range(NMAX+1):
        print(f"  b_0{k:2d} = {mp.nstr(BA[0,k], 12)}")

    print("\n--- Decay analysis ---")
    for label, B in [("Cheb-T", BT), ("Cheb-U", BU), ("Asym", BA)]:
        result = analyze_decay(B, label)
        if result:
            p_fit, data = result
            # Theoretical prediction for decay:
            # b_{0k} ≈ C · k^{-p}
            # For Cheb-T (α=β=-1/2): W^{3/2} = (1-x²)^{-3/4}
            # The integrand has endpoint exponent -3/4+1/2 = -1/4 after x=cosθ sub
            # → p = (exponent + 1) from Riemann-Lebesgue for power singularities
            print(f"    Theoretical: endpoint exponent 3α/2={3*(-0.5)/2 if 'T' in label else 3*(0.5)/2 if 'U' in label else '?'}")

    print("\n--- Generator ring analysis (Cheb-T) ---")
    find_generator(BT, BU)

    print("\nComputation complete. Tables available in:", CACHE)
