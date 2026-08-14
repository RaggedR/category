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

## ✅ Session results 2026-08-10 — sextic test DONE, conjecture graduates (refined)

Scripts + logs in repo: `jacobi_sextic.py` (exact Gram–Schmidt, 40-digit
verification of every entry against the defining integral),
`jacobi_sextic_analysis.py` (primitive-ℤ frame, integrality),
`jacobi_sextic_connection.py` (the structural theorem). α = β = −5/9.

1. **Sextic generator confirmed.** ν₀ = B(½,⅙) = 2^{2/3}πΓ(⅓)/Γ(⅔)²
   (plan's evaluation checked, 40 digits). General ratio, a=(p+2)/2:
   𝓘(p,2n)/𝓘(p,0) = (1−a)ₙ/(a)ₙ; p=−⅔ gives (⅓)ₙ/(⅔)ₙ — thirds.
2. **Two-lattice refinement of the central conjecture.** The entry
   integral (W^{3/2}) contributes the sextic/thirds lattice; the
   *polynomials themselves* (weight W) contribute a second lattice at
   eighteenths: monic Gegenbauer λ=−1/18 recurrence
   βₙ = 9n(9n−10)/((18n−1)(18n−19)) (verified n≤8), so norms live on
   ∏(9j−10), ∏(18j±1) classes. Chebyshev's single quartic lattice is the
   λ=0 degeneration (Tₙ∈ℤ[x], β=¼). Matrix = (eighteenths)ᵀ(thirds)(eighteenths).
3. **Integrality.** In the primitive-ℤ polynomial frame (PZₙ, contents =
   sliding windows of (18j−1) primes), den(Eₘₖ) = 2^a × primes from
   {5,11,17,23,…} = **6j−1 — a genuine sextic factorial family**, depending
   essentially on s=(m+k)/2 only. So an s-dependent rescale gives ℤ entries
   (same shape as Chebyshev's Q(s)); there is no clean pure-row factor.
4. **Structural theorem (the real prize; verified ALL m,k≤8 exactly).**
   With λ = α+½ (polys), μ = 3α/2+½ (i.e. W^{3/2} = (1−x²)^{μ−½}), the
   classical Gegenbauer connection formula gives every entry as a finite
   Pochhammer sum — the matching matrix is the **Gram matrix of connection
   coefficients** C^λ → C^μ in the μ-norm:
   ⟨πₘ,πₖ⟩_{W^{3/2}} = mfₘ mfₖ Σⱼ A_{(m−j)/2}(m) A_{(k−j)/2}(k) hⱼ^μ,
   A_k(n) = (λ−μ)ₖ(λ)_{n−k}(n−2k+μ)/(k!(μ)_{n−k+1}),
   hⱼ^μ/h₀^μ = (2μ)ⱼ μ/(j!(j+μ)), mfₙ = n!/(2ⁿ(λ)ₙ), and λ−μ = −α/2
   for every symmetric Jacobi weight. Row-0 closed form:
   q_{0,2n} = (2n)!(λ−μ)ₙ(λ)ₙμ/(2^{2n}(λ)_{2n} n!(μ)ₙ₊₁) — the (5/18)ₙ
   factor explains the observed primes 5,23,41,59 (≡5 mod 18).
   The Chebyshev two-term closed form ½[P(d)Q(s)/Q(d)+P(s)] should be the
   λ→0 limit. This is a *proof mechanism*: general (α,β) follows from
   Jacobi connection coefficients (also classical, single-sum).

**New pitfall (add to the list):** sympy `nsimplify` numerically re-searches
even exact Rationals and can return surd garbage (e.g. 621/14840 →
225·2^{78/145}·…). Never call nsimplify on exact quantities.

## ✅ Write-up done 2026-08-10 (evening session, commit 2ba184f)

New §7 "The Sextic Test" in `jacobi-matching.tex` (Wallis-ratio prop, ν₀
generator, two lattices, E-table + sextic denominator law, connection-Gram
theorem with DLMF §18.18(iv) cite, row-0 corollary, Chebyshev corollary);
abstract item 5, Discussion + open problem 5, code availability. Synthesis
§2.2 NOT touched — Robin's explicit choice; ask before porting anything.

**Key discovery — rationality, not limits.** The planned "λ→0 limit" does
not exist as an analytic problem: in the *monic frame* every entry
⟨πₘ,πₖ⟩/𝓘(3α+1,0) is a rational function of α (the 1/(λ)ₙ of the monic
factor cancels the (λ)_{n−ℓ} zeros in the connection coefficients).
Chebyshev is *evaluation* at α=−½; proof = two rational functions agreeing
on (−⅔,∞)∖{−½}. Verified all m,k≤8 exactly, plus the integer two-term form
and the degenerate Tₙ-in-C^{−1/4} coefficients: `jacobi_lambda0_limit.py`
(.log is gitignored by `*.log` — use `git add -f`, the repo convention).

**Honest significance assessment (agreed with Robin).** The sextic result
*deflates* the Jacobi column: connection coefficients are classical and the
generator is just the Beta-value mass, so quartic/sextic integrality is
standard special-function bookkeeping. The still-genuine mysteries are on
the Mellin side (MP Bessel-period lattice ℤA⊕ℤB, Wilson's Barnes
obstruction, Weil-form positivity) where no connection-coefficient argument
exists. Calibrate effort accordingly — Jacobi tasks below are consolidation,
not discovery.

### Next steps for the next agent (priority order)

1. **General (α,β) beta-value law — the one Jacobi task with upside.**
   W^{3/2} = (1−x)^{3α/2}(1+x)^{3β/2} is again a Jacobi weight, so the
   matching matrix is the Gram matrix of the Jacobi connection coefficients
   P^{(α,β)} → P^{(3α/2,3β/2)} (classical single sums; Askey, or
   Rainville ch. 16) in the target norm. Generator = mass
   B(3α/2+1, 3β/2+1)·2^{3(α+β)/2+1}-ish — write it once, carefully.
   The rationality trick extends: monic entries are rational functions of
   *(α,β) jointly*, so all degenerate cases (Chebyshev 3rd/4th kind,
   asymmetric ½,−½) are evaluations. Script first
   (`jacobi_general_connection.py`, exact sympy, verify against the
   already-computed asymmetric dyadic row 2C(2n,n)/−(4n+3)Cat(n)), then a
   short section or standalone note. Success criterion: a single theorem
   "entries ∈ ℚ·B(3α/2+1,3β/2+1) with denominators from explicit
   Pochhammers (−α/2)ₙ, (−β/2)ₙ, (λ)ₙ, (μ)ₙ" + the CM classification
   (which (α,β) give Γ-products collapsing to algebraic·π — the
   Wallis-power-zero family β=−α, original Task 2, becomes a corollary,
   Γ-free because B(1+t,1−t) = πt/sin(πt)·ℚ).
2. **Reverse matrices via the same mechanism (original Task 3).** Reverse
   kernel W^{1/2} ⇒ μ′ = α/2+½, λ−μ′ = +α/2 — sign-flipped Pochhammer,
   same Gram theorem. Chebyshev reverse (already in matching-coefficients
   §2.2: P̃, −Q(n−1)) should drop out as the α=−½ evaluation; check, then
   the U-reverse class-reuse question (does (−¼)ₙ/(5/4)ₙ reuse P,P̃? if
   yes, soften the "four classes exactly once" sentence in §2.2 — scope it
   to the four *forward* matrices).
3. **Ballot-path model for the Catalan row (original Task 4, the deep
   combinatorial one).** Unchanged; search literature first (Viennot/
   Flajolet weighted Motzkin paths, 4th-kind Chebyshev moments). Now also
   ask: do the *Jacobi connection coefficients* at (½,−½) have a known
   lattice-path interpretation? That would connect Tasks 1 and 4.
4. **Only if Robin asks:** port a condensed sextic paragraph to
   matching-coefficients.tex §2.2.

**Alternative worth proposing to Robin before grinding through 2–3:** pivot
to the Mellin side, where the open problems are real — e.g. prove the MP
ℤA⊕ℤB integrality (currently empirical + PSLQ) using the Bessel-moment
literature (Fresán–Sabbah–Yu), or compute the Wilson kernel's D-module rank
to test "lattice rank = ODE order". The Jacobi story is now clean enough to
leave.

## Tasks, in order (original — Task 1 done, see above)

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
