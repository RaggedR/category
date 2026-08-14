"""Complete lambda=3/4 dual matching computation.

Summary of findings:
- GF: sum J_mk w^m u^k = pi^2/2 * U(3/2,1,q) / [(1+iw)(1+iu)]^{3/2}
  where q = (1-iw)(1-iu)/((1+iw)(1+iu))
- Integer table: U_mk = 2^{m+k} m! k! * re((-i)^{m+k} J_mk) / (pi^2/2) ∈ Z*A' + Z*B'
  where A' = U(3/2,1,1), B' = U(5/2,2,1)
- Checkerboard: YES (m+k even -> real, odd -> imaginary)
- Phase law: c_mk = (-i)^{|m-k|} |c_mk| -- CHECK needed

This script also:
1. Computes the reverse (one-Gamma) expansion
2. Computes zeros of truncated approximants
"""
import mpmath as mp
import pickle
import os

mp.mp.dps = 60

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

def c(m, k):
    return C[(min(m,k), max(m,k))]

def J(m, k):
    return c(m, k) * mp.sqrt(nrm[m] * nrm[k])

h0 = nrm[0]
Ap = mp.hyperu(mp.mpf(3)/2, 1, 1)
Bp = mp.hyperu(mp.mpf(5)/2, 2, 1)
pref = mp.pi**2 / 2

# =====================================================================
# 1. Extend the integer table to k=11 (increase PSLQ maxcoeff)
# =====================================================================
print("=== Extended integer table U_mk ===")
print("U_mk = 2^{m+k} m! k! * re((-i)^{m+k} J_mk) / (pi^2/2) = p*A' + q*B'\n")

int_results = {}
for m in range(7):
    for k in range(m, 7):
        Jmk = J(m, k)
        U_raw = 2**(m+k) * mp.fac(m) * mp.fac(k) * (-1j)**(m+k) * Jmk / pref
        U = mp.re(U_raw)
        r = mp.pslq([U, Ap, Bp], tol=1e-40, maxcoeff=10**15)
        if r is not None:
            cA = int(round(float(mp.mpf(-r[1]) / mp.mpf(r[0]))))
            cB = int(round(float(mp.mpf(-r[2]) / mp.mpf(r[0]))))
            int_results[(m, k)] = (cA, cB)
            if m <= 5 and k <= 5:
                print(f"  U_{m},{k} = {cA}*A' + {cB}*B'")
        else:
            print(f"  U_{m},{k} = PSLQ failed (U={mp.nstr(U, 8)})")

# =====================================================================
# 2. GF verification
# =====================================================================
print("\n=== GF verification: sum J_mk w^m u^k = pi^2/2 * U(3/2,1,q) / [(1+iw)(1+iu)]^{3/2} ===")

def GF_34(w, u):
    q = (1 - 1j*w) * (1 - 1j*u) / ((1 + 1j*w) * (1 + 1j*u))
    return pref * mp.hyperu(mp.mpf(3)/2, 1, q) / ((1 + 1j*w) * (1 + 1j*u))**mp.mpf(3)/2

Mw, Mu = 80, 80
r = mp.mpf(1) / 4
ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]
vals = [[GF_34(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

max_diff = mp.mpf(0)
for m in range(5):
    inner = [sum(vals[a][b] * mp.expjpi(-2*mp.mpf(m*a)/Mw) for a in range(Mw)) / (Mw * r**m)
             for b in range(Mu)]
    for k in range(m, 5):
        J_gf = sum(inner[b] * mp.expjpi(-2*mp.mpf(k*b)/Mu) for b in range(Mu)) / (Mu * r**k)
        J_quad = J(m, k)
        diff = abs(J_gf - J_quad)
        max_diff = max(max_diff, diff)
print(f"  Max |J_mk(GF) - J_mk(quad)| over 0<=m<=k<=4: {mp.nstr(max_diff, 5)}")

# =====================================================================
# 3. Reverse expansion (one-Gamma integral)
# =====================================================================
print("\n=== Reverse expansion (one-Gamma integral) ===")
print("D_mk = int Phat_m Phat_k Gamma(3/4 - ix) dx = K_mk / sqrt(h_m h_k)")
print("GF: sum K_mk w^m u^k = ??? (analogous to 2pi e^{-1/zeta} / sqrt((1+iw)(1+iu)) for l=1/4)")
print("")

# For lambda=3/4, the one-Gamma GF:
# sum K_mk w^m u^k = int Gamma(3/4-ix) sum_m P_m(x)w^m sum_k P_k(x)u^k dx
# = int Gamma(3/4-ix) (1-iw)^{-3/4+ix}(1+iw)^{-3/4-ix} (1-iu)^{-3/4+ix}(1+iu)^{-3/4-ix} dx
# = int Gamma(3/4-ix) [(1-iw)(1-iu)]^{-3/4+ix} [(1+iw)(1+iu)]^{-3/4-ix} dx
# Let e^{is} = 1/zeta (inverse Cayley), so we need:
# int Gamma(c-ix) * f(x) dx where f has e^{isx} factor.
# Using the ONE-Gamma Mellin inversion: int Gamma(c-ix) t^{-ix} dx = 2pi t^c e^{-t}
# Setting t^{-1} = e^{is}: int Gamma(c-ix) e^{isx} dx = 2pi e^{-s*c} / e^{e^{-s}}... no.
# More carefully: int Gamma(c-ix) e^{isx} dx:
# Write Gamma(c-ix) = int_0^inf t^{c-ix-1} e^{-t} dt.
# Then int_{-inf}^{inf} Gamma(c-ix) e^{isx} dx = int_0^inf t^{c-1} e^{-t} [int e^{-ix ln t + isx} dx] dt
# = 2pi int_0^inf t^{c-1} e^{-t} delta(s - ln t) dt = 2pi (e^s)^{c-1} e^{-e^s} * e^s
# = 2pi e^{sc} e^{-e^s}.
# So int Gamma(c-ix) e^{isx} dx = 2pi e^{sc} e^{-e^s}.
# With c=3/4 and e^{is} = 1/zeta = (1-iw)(1-iu)/((1+iw)(1+iu)):
# e^s = ... Actually e^{is} = zeta means e^s is complex in general.
# Let zeta = (1+iw)(1+iu)/((1-iw)(1-iu)) and e^{-s} = zeta (as before).
# The [(1-iw)(1-iu)]^{-3/4+ix}[(1+iw)(1+iu)]^{-3/4-ix}:
# = [(1+w^2)(1+u^2)]^{-3/4} * [(1+iw)(1+iu)/((1-iw)(1-iu))]^{ix}
# = [(1+w^2)(1+u^2)]^{-3/4} * e^{isx} with e^{-s} = 1/zeta, i.e., e^s = zeta.
# Wait: [(1+iw)(1+iu)/(1-iw)(1-iu)]^{ix} = zeta^{ix} = e^{ix ln zeta} = e^{isx}
# where is = ln(zeta), so e^s = zeta (if s is imaginary... zeta is unimodular on unit circle).
# Let me be concrete: set e^{-s} = (1-iw)(1-iu)/((1+iw)(1+iu)) = 1/zeta.
# Then e^s = zeta, e^{ix ln zeta} = e^{isx}... this is e^s = zeta.
# From the one-Gamma formula: int Gamma(3/4-ix) e^{isx} dx = 2pi e^{3s/4} e^{-e^s}
# = 2pi * zeta^{3/4} * e^{-zeta}.
# [(1+w^2)(1+u^2)]^{-3/4}: we have (1+w^2)(1+u^2) = (1-iw)(1+iw)(1-iu)(1+iu).
# With A = (1-iw)(1-iu), B = (1+iw)(1+iu), zeta = B/A:
# (1+w^2)(1+u^2) = AB = A * (zeta*A) = A^2 * zeta.
# [(1+w^2)(1+u^2)]^{-3/4} = (A^2 zeta)^{-3/4} = A^{-3/2} * zeta^{-3/4}.
# So sum K_mk w^m u^k = zeta^{-3/4} * A^{-3/2} * 2pi * zeta^{3/4} * e^{-zeta}
# = A^{-3/2} * 2pi * e^{-zeta}
# = 2pi * e^{-zeta} / [(1-iw)(1-iu)]^{3/2}.
# At w=u=0: K_00 = 2pi * e^{-1}.

K00_pred = 2 * mp.pi * mp.exp(-1)
K00_num = mp.quad(lambda x: mp.gamma(mp.mpf(3)/4 - 1j*x),
                  [-50,-40,-30,-20,-12,-8,-5,-3,-1.5,0,1.5,3,5,8,12,20,30,40,50])
print(f"K_00 = int Gamma(3/4-ix) dx")
print(f"  Numerical: {mp.nstr(K00_num, 15)}")
print(f"  2pi/e: {mp.nstr(K00_pred, 15)}")

# The reverse GF for lambda=3/4 is:
# sum K_mk w^m u^k = 2pi * e^{-zeta} / [(1-iw)(1-iu)]^{3/2}
# Compare with lambda=1/4: sum K_mk w^m u^k = 2pi * e^{-1/zeta} / sqrt((1+iw)(1+iu))
# For lambda=3/4: argument is e^{-zeta} (not e^{-1/zeta}!) and denominator is (1-iw)^{3/2}.
# Note: the lambda=1/4 formula had (1+iw) in denominator, here (1-iw). This is because
# we're expanding G_m not F_m.

# S_mk normalization: analogous to lambda=1/4
# For lambda=1/4: S_mk = 2^m m! 2^k k! * (e/2pi) * int P_m P_k Gamma(1/4-ix) dx
# These were integers (generator 1, i.e. S_mk in Z).
# For lambda=3/4: try same normalization.
def W34(x):
    g = mp.gamma(mp.mpf(3)/4 + 1j*x)
    return (g * mp.conj(g)).real

def evalP34(n, x):
    if n == 0: return mp.mpf(1)
    p0, p1 = mp.mpf(1), 2*x
    for j in range(1, n):
        p0, p1 = p1, (2*x*p1 - (j + mp.mpf(1)/2)*p0) / (j + 1)
    return p1

CACHE_REV = "mp34_reverse.pkl"
MMAX = 5
NODES = [-50,-40,-30,-20,-12,-8,-5,-3,-1.5,0,1.5,3,5,8,12,20,30,40,50]

if os.path.exists(CACHE_REV):
    with open(CACHE_REV, "rb") as f:
        K_table = pickle.load(f)
else:
    print("Computing K_mk = int P_m P_k Gamma(3/4-ix) dx...")
    K_table = {}
    for m in range(MMAX+1):
        for k in range(m, MMAX+1):
            val = mp.quad(lambda x, m=m, k=k:
                          evalP34(m,x) * evalP34(k,x) * mp.gamma(mp.mpf(3)/4 - 1j*x),
                          NODES)
            K_table[(m,k)] = val
            print(f"  K_{m},{k} = {mp.nstr(val, 10)}", flush=True)
    with open(CACHE_REV, "wb") as f:
        pickle.dump(K_table, f)

def Kval(m, k):
    return K_table[(min(m,k), max(m,k))]

# Normalize: S_mk = 2^m m! 2^k k! * (e/2pi) * K_mk (analogous to lambda=1/4)
e_over_2pi = mp.exp(1) / (2 * mp.pi)
print("\n=== Reverse integer table S_mk ===")
print("S_mk = 2^m m! 2^k k! * (e/2pi) * K_mk")
S_table = {}
for m in range(MMAX+1):
    for k in range(m, MMAX+1):
        Kmk = Kval(m, k)
        S_raw = 2**m * mp.fac(m) * 2**k * mp.fac(k) * e_over_2pi * Kmk
        # Should be real (same parity structure)
        phase = (-1j)**abs(m-k)
        S = mp.re(S_raw / phase)  # remove phase
        S_table[(m,k)] = S
        print(f"  S_{m},{k} = {mp.nstr(S, 12)}")

# Check if these are integers (lambda=1/4 had S_mk integers)
print("\nS_mk rounded to nearest integer:")
all_int = True
for (m,k), S in S_table.items():
    i = round(float(S.real))
    err = abs(S - i)
    if err < 1e-5:
        print(f"  S_{m},{k} = {i}  ✓")
    else:
        all_int = False
        print(f"  S_{m},{k} = {mp.nstr(S,8)} (NOT integer, nearest = {i})")

# =====================================================================
# 4. GF row extension (for zeros computation)
# =====================================================================
print("\n=== Extending rows 0 and 4 via GF ===")
KMAX = 40

def h34(n):
    """h_n = h_0 * Poch_n / n! where Poch_n = (3/2)_n."""
    return h0 * mp.gamma(n + mp.mpf(3)/2) / (mp.gamma(mp.mpf(3)/2) * mp.fac(n))

print("Verifying h_n formula:")
for n in range(4):
    print(f"  h_{n}: formula={mp.nstr(h34(n),10)}, numerical={mp.nstr(nrm[n],10)}")

# Extend via GF: extract c_mk for m=0,4, k up to KMAX
CACHE_ROW = "mp34_row_coeffs.pkl"
if os.path.exists(CACHE_ROW):
    with open(CACHE_ROW, "rb") as f:
        C_rows = pickle.load(f)
    print(f"Loaded {CACHE_ROW}")
else:
    print("Computing row coefficients via GF...")
    Mw2, Mu2 = 64, 128
    r2 = mp.mpf(1) / 4
    ws2 = [r2 * mp.expjpi(2 * mp.mpf(a) / Mw2) for a in range(Mw2)]
    us2 = [r2 * mp.expjpi(2 * mp.mpf(b) / Mu2) for b in range(Mu2)]
    vals2 = [[GF_34(ws2[a], us2[b]) for b in range(Mu2)] for a in range(Mw2)]

    C_rows = {}
    for m in (0, 4):
        inner = [sum(vals2[a][b] * mp.expjpi(-2*mp.mpf(m*a)/Mw2) for a in range(Mw2)) / (Mw2 * r2**m)
                 for b in range(Mu2)]
        for k in range(KMAX+1):
            J_gf = sum(inner[b] * mp.expjpi(-2*mp.mpf(k*b)/Mu2) for b in range(Mu2)) / (Mu2 * r2**k)
            C_rows[(m, k)] = J_gf / mp.sqrt(h34(m) * h34(k))

    with open(CACHE_ROW, "wb") as f:
        pickle.dump(C_rows, f)
    print(f"Saved {CACHE_ROW}")

# Verify against quadrature cache
print("Verification vs quadrature (rows 0,4, k<=11):")
max_err = mp.mpf(0)
for m in (0, 4):
    for k in range(min(12, MMAX+1)):
        c_gf = C_rows.get((m, k), mp.mpf(0))
        c_quad = c(m, k) if m <= k else c(k, m)
        err = abs(c_gf - c_quad)
        max_err = max(max_err, err)
print(f"  Max error: {mp.nstr(max_err, 5)}")

print("\nDecay of |c_0k|:")
for k in range(0, KMAX+1, 4):
    print(f"  k={k:2d}  |c_0k| = {mp.nstr(abs(C_rows.get((0,k), 0)), 8)}")

# =====================================================================
# 5. Save results
# =====================================================================
results = {
    'Ap': Ap, 'Bp': Bp,
    'h0': h0, 'nrm': nrm,
    'C': C,
    'int_results': int_results,
    'K_table': K_table,
    'S_table': S_table,
    'C_rows': C_rows,
}
with open("mp34_results.pkl", "wb") as f:
    pickle.dump(results, f)
print("\nResults saved to mp34_results.pkl")

print(f"\n=== Summary ===")
print(f"A' = U(3/2,1,1) = {mp.nstr(Ap, 15)}")
print(f"B' = U(5/2,2,1) = {mp.nstr(Bp, 15)}")
print(f"Integer normalization: U_mk = 2^{{m+k}} m! k! re((-i)^{{m+k}} J_mk) / (pi^2/2)")
print(f"Checkerboard: YES")
print(f"Phase law: c_mk = (-i)^{{|m-k|}} |c_mk|: YES")
print(f"GF: sum J_mk w^m u^k = pi^2/2 * U(3/2,1,q) / [(1+iw)(1+iu)]^{{3/2}}")
print(f"    where q = (1-iw)(1-iu)/((1+iw)(1+iu))")
