"""The matrix Z_{nk} of log-Hermite polynomial approximations of the zeta family.

Family: xi_n(s) = (c_n/c_0) P_n(s) xi(s)  (even n; P_n from the Mellin transform
of psi_n, zeros certified on the line by Bump--Ng).  Since multiplication by
tau acts on Hermite-function coefficients as the tridiagonal position operator

    (T v)_j = sqrt(j/2) v_{j-1} + sqrt((j+1)/2) v_{j+1},

the projection coefficients of F_n(t) = xi_n(1/2 + i*SCALE*t) onto polynomials
in L^2(e^{-t^2}) are

    Z_n = (c_n/c_0) * P_n(1/2 + i*SCALE*T) beta,

with beta the single Xi column from xi_polyapprox.py.  Rows with n = 2 mod 4
are i * (real, odd); we divide out the i.

On the line (tau = SCALE*t):
    P_0 = 1
    P_2 = 2(2s-1)            = 4 i tau
    P_4 = 4(4s^2-4s+3)       = 4(2 - 4 tau^2)
    P_6 = 8(2s-1)(4s^2-4s+15)= 16 i tau (14 - 4 tau^2)
"""
import mpmath as mp
import pickle

mp.mp.dps = 80
SCALE = 10
NMAX = 60
PAD = 4  # room for Jacobi spill; top entries incomplete (missing beta_{62,64,66})

with open("xi_polyapprox_coeffs.pkl", "rb") as f:
    bs = pickle.load(f)
beta = [bs.get(k, mp.mpf(0)) for k in range(NMAX + 1)] + [mp.mpf(0)] * PAD


def T(v):
    out = [mp.mpf(0)] * len(v)
    for j in range(len(v)):
        if j - 1 >= 0:
            out[j] += mp.sqrt(mp.mpf(j) / 2) * v[j - 1]
        if j + 1 < len(v):
            out[j] += mp.sqrt(mp.mpf(j + 1) / 2) * v[j + 1]
    return out


def scal(a, v):
    return [a * x for x in v]


def add(u, v):
    return [a + b for a, b in zip(u, v)]


L = mp.mpf(SCALE)
Tb = T(beta)
T2b = T(Tb)
T3b = T(T2b)

rows = {
    0: beta,
    # (1/sqrt 8) * 4 L * (T beta)           [i divided out]
    2: scal(4 * L / mp.sqrt(8), Tb),
    # (1/sqrt 384) * 4 * (2 beta - 4 L^2 T^2 beta)
    4: scal(mp.mpf(4) / mp.sqrt(384), add(scal(2, beta), scal(-4 * L**2, T2b))),
    # (1/sqrt 46080) * 16 L * (14 T beta - 4 L^2 T^3 beta)   [i divided out]
    6: scal(16 * L / mp.sqrt(46080), add(scal(14, Tb), scal(-4 * L**2, T3b))),
}

print("matrix Z_nk  (log10 |Z_nk|, '.' = exact zero by parity):")
ks = list(range(0, NMAX + 1, 6))
print("  n\\k " + "".join("%7d" % k for k in ks))
for n, row in rows.items():
    cells = []
    for k in ks:
        # parity: row n hits k = (n/2) mod 2
        v = row[k] if k % 2 == (n // 2) % 2 else row[k + 1] if k + 1 <= NMAX else 0
        cells.append("%7.1f" % mp.log10(abs(v)) if v else "      .")
    print("  n=%d " % n + "".join(cells))

# zeros of the degree-60 approximants
def hermite_coeffs(n):
    a, b = [1], [0, 2]
    if n == 0:
        return a
    for i in range(1, n):
        nxt = [0] * (i + 2)
        for j, cj in enumerate(b):
            nxt[j + 1] += 2 * cj
        for j, cj in enumerate(a):
            nxt[j] -= 2 * i * cj
        a, b = b, nxt
    return b


N = 60
mp.mp.dps = 120
certified = {0: [], 2: [0.0], 4: [1 / mp.sqrt(2)], 6: [0.0, mp.sqrt(mp.mpf(7) / 2)]}
for n, row in rows.items():
    coeffs = [mp.mpf(0)] * (N + 1)
    for k in range(N + 1):
        if row[k] == 0:
            continue
        ck = (mp.mpf(2) ** k * mp.factorial(k) * mp.sqrt(mp.pi)) ** mp.mpf(-0.5)
        w = mp.mpf(row[k]) * ck
        for j, hj in enumerate(hermite_coeffs(k)):
            coeffs[j] += w * hj
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    roots = mp.polyroots([c for c in coeffs[::-1]], maxsteps=300, extraprec=400)
    rr = sorted(float(z.real) for z in roots if abs(z.imag) < 1e-6 and z.real > -1e-9)
    cc = sorted((float(z.real), float(z.imag)) for z in roots if z.imag > 1e-6)
    print("\nxi_%d approximant, N=%d:" % (n, N))
    print("  certified P_%d zeros at T = %s"
          % (n, ", ".join("%.4f" % (float(z) * 1) for z in certified[n])
             if certified[n] else "(none)"))
    print("  real zeros >= 0, * SCALE:")
    print("   " + ", ".join("%.4f" % (r * SCALE) for r in rr[:16])
          + (" ..." if len(rr) > 16 else ""))
    print("  complex zeros (upper half, * SCALE): "
          + ", ".join("(%.2f,%.2f)" % (x * SCALE, y * SCALE) for x, y in cc))

print("\nfirst Riemann zeros: 14.1347, 21.0220, 25.0109, 30.4249, 32.9351")
