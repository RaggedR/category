# Plan: General Jacobi Matching Matrices

*Handover for the next agent. Written 2026-08-10 after commits e95fe23 (phase-error fix) and 04f2812 (U + asymmetric integrality).*

First: read `~/.claude/AGENT.md`, then the memory index at
`~/.claude/projects/-Users-robin-git-category/memory/MEMORY.md`
(especially `project_matching_synthesis.md`). The papers are
`matching-coefficients.tex` (THE synthesis paper, §2.2 Chebyshev) and
`jacobi-matching.tex` (standalone, corrected).

## The construction (one sentence each)

Kernel W → orthonormal polynomials P̂ₙ → functions Φₙ = P̂ₙ√W.
Forward matching coefficients bₘₖ = ∫P̂ₘP̂ₖW^{3/2}dx (these *approximate*
Φₘ by polynomials — say "approximation", never "expansion", for the forward
direction). Reverse coefficients use W^{1/2}.

## Established results (all quadrature-verified, do not re-derive from scratch)

Quartic factorials: P(n)=1·5···(4n−3), Q(n)=3·7···(4n−1),
P̃(n)=5·9···(4n+1), Q⁺(n)=7·11···(4n+3). Write s=(m+k)/2, d=|m−k|/2.

| matrix | kernel | rescale by | entry | edge | diagonal |
|---|---|---|---|---|---|
| forward T | (1−x²)^{−1/2} | Q(s)/J(0) | ½[P(d)Q(s)/Q(d)+P(s)] ∈ ℤ₊ | P(n) = A007696 | ½[Q+P] = 1,2,13,138 ∉ OEIS |
| reverse T | same, W^{1/2} | P̃(s)/J̃(0) | ∈ ℤ, signs alternate | −Q(n−1) | ½[P̃−Q(n−1)] = 1,2,21,282 ∉ OEIS |
| forward U | (1−x²)^{+1/2} | 3Q⁺(s+1)/J(0) | ½[Q⁺(s+1)ρ(d)+3P(s)] ∈ ℤ, ρ(n)=−3P(n−1)/Q⁺(n) | −15P(n−1) | ½[Q⁺(n+1)+3P(n)] = 5,40,585 ∉ OEIS |
| asymmetric (½,−½) | (1−x)^{1/2}(1+x)^{−1/2} | 2^k√2 | ℤ[½]·√2, whole matrix dyadic | row 0: 2C(2n,n) even, −(4n+3)Cat(n) odd | — |

Generator constant J(0) = √π Γ(¼)/Γ(¾) = Γ(¼)²/√(2π) — the lemniscatic
period (CM curve y² = 1−x⁴). All odd-parity entries of symmetric kernels
vanish. Odd entries of the asymmetric kernel do NOT vanish.

Key mechanism: x=cosθ collapses TₘTₖ to two cosines; the Wallis integral is
𝓘(p,q) = ∫₀^π sinᵖθ cos(qθ)dθ = **cos(πq/2)** · πΓ(p+1) / (2ᵖ Γ((p+q+2)/2) Γ((p−q+2)/2)).

## ⚠ Pitfalls (each one burned a previous session)

1. **The phase trap.** The cos(πq/2) factor in the Wallis formula is
   essential. Dropping it creates phantom odd-parity entries and sign errors
   that *look like structure* (a fake second generator γ_o lived in
   jacobi-matching.tex for days). Always verify closed forms against
   `mp.quad` of the *defining integral*, never against a re-derivation of the
   same combinatorics.
2. **Quadrature at endpoint singularities.** x-space quadrature of
   (1−x²)^{−3/4} is unreliable near ±1. Substitute θ: integrate
   `mp.quad(lambda t: ..., [0, mp.pi])` at `mp.dps = 40`, splitting at
   π/4, π/2, 3π/4 for stubborn cases.
3. **Rational recognition.** Use `Fraction(float(v)).limit_denominator(10**6)`
   and then assert agreement to ≥ 10 digits; for two-generator lattices use
   mpmath PSLQ with `maxcoeff ≥ 1e10`.
4. **Robin's conventions.** Real theorem environments with his statements
   verbatim; "triangular" matrix (not "symmetric") if he says so; forward =
   approximation, reverse = expansion; no colons mid-sentence in prose;
   unicode math in conversation, LaTeX only in .tex.

## The central conjecture (why this matters)

**Matching arithmetic = period arithmetic.** The generator ring of the
matching matrix appears to be the period lattice of the weight:
Gaussian → ℤ; Chebyshev → quartic Pochhammers (¼)ₙ, (¾)ₙ over the
lemniscatic period Γ(¼)²; asymmetric half-integer → dyadic Catalan/central
binomial (ballot paths); MP → Bessel-period lattice ℤA⊕ℤB.

## Tasks, in order

### 1. The sextic test (decisive, do first)

Prediction: a symmetric Jacobi weight whose Wallis power p = 3α+1 has
(p+2)/2 ∈ ⅓+ℤ should produce **sextic** structure — Γ(⅓)-type generator
(period of y² = 1−x⁶, CM by ℤ[ω]) and Pochhammer ratios stepping by
thirds. Candidate: α = β = −5/9 gives p = −⅔, so
𝓘(p,0) = πΓ(⅓)/(2^{−2/3}Γ(⅔)²) — check this evaluation yourself.
- The polynomials are Gegenbauer C^{(λ)}, λ = α+½ = −1/18; they do NOT
  collapse to single cosines. Expand them in the Chebyshev basis (exact
  rational coefficients via recurrence) so every entry becomes a finite
  ℚ-combination of 𝓘(−⅔, 2n).
- First derive 𝓘(p,2n)/𝓘(p,0) as a ratio of shifted Pochhammers for
  general p (reflection formula; watch the phase-sign cancellation).
- Compute rows 0–4 exactly. Question: is there a single row factor making
  entries integral? Does a "sextic factorial" family appear?
- If yes: the conjecture graduates. If no: report *how* it fails — a
  rational-but-not-integral answer is still informative.

### 2. The Wallis-power-zero family (cheap, surprising if true)

(α,β) = (a,−a) always has p = 3(a−a)/2·… — check: does every member of
this one-parameter family have Γ-free (rational or algebraic) entries, or
is a = ½ special because the polynomials are trigonometric (3rd/4th-kind
Chebyshev)? Test a = ¼ numerically. Also do the mirror case
(−½, ½) (third-kind Vₙ, cos((n+½)θ)/cos(θ/2)) — expect the same
sequences up to signs; confirm.

### 3. Reverse matrices for U and asymmetric (bookkeeping + one honest check)

Reverse-U uses W^{1/2} = (1−x²)^{1/4}, effective p = ½, whose Pochhammer
ratios are (−¼)ₙ/(5/4)ₙ — this REUSES classes P and P̃. That would break
the "four classes each exactly once" sentence at the end of §2.2 of
matching-coefficients.tex. Compute it, and if the reuse is real, soften
that sentence (it currently refers only to the four matrices computed
so far — make the scope explicit or restate the poetry honestly).
Reverse asymmetric: expect dyadic again; find the row-0 combinatorics.

### 4. A ballot-path model for the Catalan row (the deep one)

Row 0 of the asymmetric kernel is 2C(2n,n) / −(4n+3)Cat(n). Catalan
numbers demand a bijection. Find a weighted-lattice-path or matching model
whose transfer matrix reproduces the full dyadic matrix, in the style of
§4 (Automata) of the synthesis paper — states, transitions, weights,
start/end rule. Search first (Robin's rule): fourth-kind Chebyshev moment
combinatorics, Viennot/Flajolet weighted Motzkin paths for W = tan(θ/2)
weights, "central binomial + Catalan interleaved" in OEIS/literature.

### 5. Only then: general (α,β)

With data points at exponents {0, ±½, −5/9}, guess the general law for
the generator ring as a function of (α,β) (expect: generated by the beta
value B(3α/2+…, 3β/2+…), i.e. a period of y^N = xᵃ(1−x)ᵇ). State it as a
conjecture in the synthesis paper only when at least three independent
exponents confirm.

## Success criteria

- Every claimed closed form checked against 40-digit quadrature of the
  defining integral for m,k ≤ 8 before it enters any .tex file.
- New sequences checked against OEIS (JSON API: oeis.org/search?q=...&fmt=json;
  a bare `null` response means no match).
- tectonic compiles clean; commit to ~/git/category (repo exists, main branch).
- Update `project_matching_synthesis.md` memory + MEMORY.md index.
- If an approach fails three times, stop and write down why (three-strike rule).
