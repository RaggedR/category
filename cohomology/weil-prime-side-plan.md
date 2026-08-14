# Plan: The Prime Side in the MP Basis

**Goal.** Theorem 5.1 of `matching-coefficients.tex` shows the F-basis diagonalizes the
archimedean term of the Weil explicit formula, with an integer lattice ℤA⊕ℤB on the
Γ-side. Test the converse: push the prime side −2ΣΛ(n)n^{−½}g(log n) through the same
basis and see whether the primes compress the way the Γ-side did. If they do, the paper's
infrastructure claim upgrades to a structural one.

**Starting point.** `weil_form_mp.py` — explicit formula validated to 21 digits against
the first 60 zeros (stage 1); basis kernel h_mn = (r²+¼)P̂ₘP̂ₙ|Γ(¼+ir/2)|², x = r/2.
Pole bug already fixed: (r²+¼) does NOT kill h(±i/2); residue = 2Γ(½)P̂ₘP̂ₙ(±i/4).

## Stage 0 — Closed form for g_mn (the real work)

The prime side needs g_mn(u) = (1/2π)∫h_mn(r)e^{−iur}dr. Do this **exactly**, not
numerically:

- Known: ∫_ℝ |Γ(λ+ix)|² e^{−2iux} dx = π Γ(2λ) / (cosh u)^{2λ} (verify constant for our
  x = r/2 convention). For λ = ¼ the kernel is (cosh u)^{−1/2}.
- P̂ₘP̂ₙ(r/2) is a polynomial in r → acts as a differential operator in u on that kernel;
  (r²+¼) likewise (−d²/du² + ¼).
- Derivatives of (cosh u)^{−1/2} stay in ℚ[tanh u, sech u]·(cosh u)^{−1/2}.
- At u = log n: cosh(log n) = (n+1/n)/2, tanh(log n) = (n²−1)/(n²+1) — **algebraic in n**.

Expected output: P_mn = 2 Σ Λ(n) n^{−½} · R_mn(n) with R_mn ∈ ℚ(√n, √(n²+1))-type
explicit algebraic functions. Verify each g_mn(u) against numerical quadrature (10+ digits,
several u) before trusting.

## Stage 1 — Compute P_mn to high precision

Using the closed form (fast, no oscillatory quadrature): dps = 50, m,n ≤ 6, both parity
blocks. Convergence: tail beyond cutoff N is ~ Σ_{n>N} Λ(n) n^{−½} R_mn(n); R_mn decays
like n^{−¼}·(polynomially bounded) — estimate the tail honestly, push cutoff until 40
stable digits. If the closed form gives per-prime terms, also try acceleration via
ψ(x) = x + error (replace Λ-sum by its smooth part + zeros — but that reintroduces zeros;
keep as a check only).

## Stage 2 — Structure hunt on P_mn

In order of strength:
1. **PSLQ** P_mn against {A, B}×{1, log2, logπ, π, γ} and against each other
   (row/column ℚ-relations). Lattice constants: A,B = √(2e)·K₀,K₁(½).
2. Sign pattern / monotonicity / bandedness of the matrix.
3. Compare magnitude with arch side A_mn: is W = pole + arch − prime positive
   definite as m grows (extend stage-2 blocks of the script beyond m ≤ 4)?
4. Per-prime structure: are the R_mn(n) values themselves in a lattice over n
   (e.g. rows of R_mn as polynomials in 1/(n+1))?

## Success / failure criteria

- **Compression** (relations found at 2× dps confirmation): new section for the paper;
  the basis sees arithmetic.
- **No compression**: still publishable as a remark — primes are generic in the F-basis,
  the integrality is purely archimedean (consistent with the Hahn/Wilson D-module
  reading). Either outcome is informative; say which one honestly.

## Pitfalls (learned the hard way — see memory)

- **PSLQ**: maxcoeff ≥ 1e10 always (Hahn lesson); any relation found at dps d must be
  re-verified at 2d or it is untrusted even with small coefficients (Wilson lesson).
- Working precision vs number of PSLQ inputs: ~12 digits per input value is the
  spurious-relation threshold.
- Parity: h_mn is even iff m ≡ n (mod 2); odd-parity entries pair with nothing — the
  form is block-diagonal, don't compute cross-parity prime terms.
- MP norm convention: h_m = 2πΓ(m+½)/(√2·m!) (φ = π/2, over ℝ) — match the script.
- Robin reads unicode math in conversation, LaTeX only in .tex files.

First, read ~/.claude/AGENT.md if you are a sub-agent. Main references:
`matching-coefficients.tex` §5, `weil_form_mp.py`, memory file
`project_matching_synthesis.md`.
