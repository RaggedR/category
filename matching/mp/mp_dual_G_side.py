"""G-side coefficient tables for the two Meixner--Pollaczek dual pairs.

The paper (matching-coefficients.tex) tabulates, for lambda = 1/4 and 3/4,
the F-side matrices only:

  forward  c+_mk = int Phat_m Phat_k Gamma(l+ix)^2 Gamma(l-ix) dx
  reverse  D_mk  = int Phat_m Phat_k Gamma(l-ix)   dx

Here we compute the G-side analogues DIRECTLY (no conjugation shortcut):

  forward  c-_mk = int Phat_m Phat_k Gamma(l+ix) Gamma(l-ix)^2 dx
  reverse  Dt_mk = int Phat_m Phat_k Gamma(l+ix) dx

rescale them exactly as the paper rescales the F-side, PSLQ the forward
entries against (A, B) resp. (A', B'), round the reverse entries to
integers, and compare with the published F-side tables.
"""
import mpmath as mp

mp.mp.dps = 40
NODES = [-60, -50, -40, -30, -20, -12, -8, -5, -3, -1.5, 0,
         1.5, 3, 5, 8, 12, 20, 30, 40, 50, 60]

I = mp.mpc(0, 1)


def evalP(lam, n, x):
    """MP polynomial P_n^{(lam)}(x; pi/2):
    (j+1) P_{j+1} = 2x P_j - (j - 1 + 2*lam) P_{j-1}, P_0=1, P_1=2*lam...
    NB the paper's normalisation has P_1 = 2x for both families and
    recurrence coefficient (j - 1/2) for lam=1/4, (j + 1/2) for lam=3/4,
    i.e. (j - 1 + 2*lam)."""
    if n == 0:
        return mp.mpf(1)
    p0, p1 = mp.mpf(1), 2 * x
    for j in range(1, n):
        p0, p1 = p1, (2 * x * p1 - (j - 1 + 2 * lam) * p0) / (j + 1)
    return p1


def quadP(lam, m, k, kern):
    return mp.quad(lambda x: evalP(lam, m, x) * evalP(lam, k, x) * kern(x),
                   NODES)


# ---------------- published F-side tables ----------------
T14_F = {(0, 0): (1, 0), (0, 1): (3, -2), (0, 2): (-23, 12),
         (0, 3): (-241, 138), (0, 4): (3265, -1800), (0, 5): (52907, -29706),
         (1, 1): (-21, 12), (1, 2): (-205, 114), (1, 3): (2575, -1440),
         (1, 4): (39411, -21978), (1, 5): (-707389, 394932),
         (2, 2): (2369, -1320), (2, 3): (33743, -18834),
         (2, 4): (-569719, 317892), (2, 5): (-11092229, 6190602),
         (3, 3): (-529549, 295572), (3, 4): (-9675553, 5399634)}

S14_F = {(0, 0): 1, (0, 1): 3, (0, 2): 5, (0, 3): 11, (0, 4): 249,
         (0, 5): -2621, (1, 1): 7, (1, 2): 47, (1, 3): 399, (1, 4): -2005,
         (1, 5): 263, (2, 2): 473, (2, 3): -617, (2, 4): 22109,
         (2, 5): -402065, (3, 3): 31511, (3, 4): -412557, (3, 5): -3017969,
         (4, 4): -1916623, (4, 5): -29235989, (5, 5): -1012892921}

U34_F = {(0, 0): (1, 0), (0, 1): (3, -6), (0, 2): (39, -60),
         (0, 3): (513, -870), (0, 4): (8673, -14280), (1, 1): (33, -60),
         (1, 2): (453, -750), (1, 3): (7035, -11760),
         (1, 4): (128643, -213630), (2, 2): (6609, -10920),
         (2, 3): (112695, -187350), (2, 4): (2215047, -3679860),
         (3, 3): (2078169, -3456420), (3, 4): (43522881, -72346830)}

S34_F = {(0, 0): 1, (0, 1): -1, (0, 2): -9, (0, 3): -43, (0, 4): -95,
         (0, 5): 1671, (1, 1): 15, (1, 2): 23, (1, 3): -283, (1, 4): -4767,
         (1, 5): -52649, (2, 2): 529, (2, 3): 5475, (2, 4): 31703,
         (2, 5): -65119, (3, 3): -15145, (3, 4): 363979, (3, 5): 8950653,
         (4, 4): -8709567, (4, 5): -100346585, (5, 5): 123082079}


def phase_split(z):
    """Return (phase, real value) with phase in {1, i, -1, -i} chosen so
    that z = phase * value with value real (to working precision)."""
    for ph in (mp.mpc(1), I, mp.mpc(-1), -I):
        v = z / ph
        if abs(mp.im(v)) < mp.mpf(10) ** (-20) * max(1, abs(v)):
            return ph, mp.re(v)
    raise ValueError(f"no clean phase for {z}")


def report_forward(lam, h0, gens, gen_names, F_table, kmax=5, mmax=3):
    a, b = gens
    print(f"\n=== lambda = {lam}: FORWARD G-side "
          f"T-_mk = 2^m m! 2^k k! / h0 * int P P G(l+ix) G(l-ix)^2 ===")
    for m in range(mmax + 1):
        for k in range(m, kmax + 1):
            if (m, k) not in F_table:
                continue
            Jm = quadP(lam, m, k, lambda x: mp.gamma(lam + I * x)
                       * mp.gamma(lam - I * x) ** 2)
            T = (2 ** m * mp.fac(m) * 2 ** k * mp.fac(k)) / h0 * Jm
            ph, Tre = phase_split(T)
            rel = mp.pslq([Tre, a, b], tol=mp.mpf(10) ** (-25),
                          maxcoeff=10 ** 9)
            ca, cb = F_table[(m, k)]
            resid = abs(Tre - (ca * a + cb * b))
            resid_m = abs(Tre + (ca * a + cb * b))
            match = ("= F entry" if resid < 1e-18 else
                     "= -(F entry)" if resid_m < 1e-18 else "DIFFERENT")
            phs = {1: "1", -1: "-1"}.get(int(mp.re(ph)) if ph.imag == 0
                                         else None,
                                         "i" if ph == I else "-i")
            print(f"  ({m},{k}) phase {phs:>2}: PSLQ {rel}   "
                  f"F-side ({ca}{gen_names[0]}+{cb}{gen_names[1]})  -> {match}")


def report_reverse(lam, F_table, nmax=5):
    print(f"\n=== lambda = {lam}: REVERSE G-side "
          f"St_mk = 2^m m! 2^k k! (e/2pi) int P P Gamma(l+ix) dx ===")
    for m in range(nmax + 1):
        for k in range(m, nmax + 1):
            Dm = quadP(lam, m, k, lambda x: mp.gamma(lam + I * x))
            S = (2 ** m * mp.fac(m) * 2 ** k * mp.fac(k)) * mp.e \
                / (2 * mp.pi) * Dm
            ph, Sre = phase_split(S)
            n = mp.nint(Sre)
            err = abs(Sre - n)
            phs = "1" if ph == 1 else "-1" if ph == -1 else \
                  "i" if ph == I else "-i"
            Sf = F_table[(m, k)]
            rel = "= F" if int(n) == Sf else \
                  "= -F" if int(n) == -Sf else "DIFFERENT"
            print(f"  ({m},{k}) phase {phs:>2}: {int(n):>12d}  "
                  f"(int err {mp.nstr(err, 3)})  F-side {Sf:>12d} -> {rel}")


def main():
    lam14, lam34 = mp.mpf(1) / 4, mp.mpf(3) / 4

    # generators
    A = mp.sqrt(2 * mp.e) * mp.besselk(0, mp.mpf(1) / 2)
    B = mp.sqrt(2 * mp.e) * mp.besselk(1, mp.mpf(1) / 2)
    Ap = mp.hyperu(mp.mpf(3) / 2, 1, 1)
    Bp = mp.hyperu(mp.mpf(5) / 2, 2, 1)

    h0_14 = mp.sqrt(2) * mp.pi ** mp.mpf(1.5)     # int |G(1/4+ix)|^2 dx
    # lambda=3/4 rescaling divides by pi^2/2 (paper's U_mk), not by h0
    h0_34 = mp.pi ** 2 / 2

    # sanity: J- = conj(J+) spot check
    Jp = quadP(lam14, 0, 1, lambda x: mp.gamma(lam14 + I * x) ** 2
               * mp.gamma(lam14 - I * x))
    Jm = quadP(lam14, 0, 1, lambda x: mp.gamma(lam14 + I * x)
               * mp.gamma(lam14 - I * x) ** 2)
    print(f"conjugation check |J-(0,1) - conj(J+(0,1))| = "
          f"{mp.nstr(abs(Jm - mp.conj(Jp)), 3)}")

    report_forward(lam14, h0_14, (A, B), ("A", "B"), T14_F)
    report_forward(lam34, h0_34, (Ap, Bp), ("A'", "B'"), U34_F, kmax=4)
    report_reverse(lam14, S14_F)
    report_reverse(lam34, S34_F)


if __name__ == "__main__":
    main()
