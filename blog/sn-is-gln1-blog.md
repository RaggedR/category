# I Only Know Two Finite Groups

I only know two finite groups: S_n and GL_n(F_q).

Actually I only *know* one finite group: S_n. Because all sorts of number theory sneaks into GL_n(F_q) and I'm scared of number theory.

This is a confession, but it's also a research programme. Can we think of S_n as GL_n(F_1) — the general linear group over the "field with one element"?

There is no field with one element. A field needs 0 and 1, and they must be distinct. One element isn't enough. But the answer to the question is yes, and the reason is one of the most beautiful accidents in mathematics.

## The two groups

The symmetric group S_n is the group of permutations of {1, 2, ..., n}. It has n! elements. Every finite group embeds into some S_n (Cayley's theorem), so in a sense S_n is the universal finite group. It's also the group you meet first, the one you can draw pictures of, the one whose representation theory is governed by Young tableaux and doesn't require you to think about roots of unity.

The general linear group GL_n(F_q) is the group of invertible n×n matrices over the finite field with q elements, where q = p^k for some prime p. Its order is

|GL_n(F_q)| = (q^n − 1)(q^n − q)(q^n − q²)⋯(q^n − q^{n−1})

which counts: choose the first row (any nonzero vector: q^n − 1 options), then the second (anything not in the span of the first: q^n − q options), and so on. You're building a basis one vector at a time, and each step removes a dimension's worth of choices.

This is a beautiful formula. It's also where the number theory starts creeping in. The character theory of GL_n(F_q) involves Deligne–Lusztig theory, étale cohomology, Frobenius endomorphisms. The representation ring depends on the arithmetic of F_q. I understand the formula for the order. I do not understand the representations. This is what I mean by scared.

S_n is combinatorics. GL_n(F_q) is arithmetic geometry in disguise. But they're not unrelated.

## The first hint: the Weyl group

GL_n(F_q) contains a copy of S_n. The permutation matrices — the ones with exactly one 1 in each row and column, everything else 0 — form a subgroup isomorphic to S_n. This is the **Weyl group**.

More than that: GL_n(F_q) is *built from* S_n. The Bruhat decomposition says

GL_n(F_q) = ∐_{w ∈ S_n} BwB

where B is the Borel subgroup (upper triangular matrices). The big group is stratified by the small group. Every matrix in GL_n(F_q) is determined by a permutation w ∈ S_n (telling you the combinatorial shape) plus some data from B (telling you the arithmetic). S_n is the skeleton. The field dresses it.

So S_n isn't just sitting inside GL_n(F_q). It's the *structural core*. When you remove the arithmetic, the permutation is what's left.

## The q → 1 miracle

Now look at what happens to the order formula as q → 1.

|GL_n(F_q)| = (q^n − 1)(q^n − q)(q^n − q²)⋯(q^n − q^{n−1})

You can factor out: pull q^{0+1+2+⋯+(n−1)} = q^{n(n−1)/2} from the product, giving

|GL_n(F_q)| = q^{n(n−1)/2} · (q^n − 1)(q^{n−1} − 1)⋯(q − 1)

The second factor involves q-integers. Define [k]_q = (q^k − 1)/(q − 1) = 1 + q + q² + ⋯ + q^{k−1}. Then

(q^n − 1)(q^{n−1} − 1)⋯(q − 1) = (q − 1)^n · [n]_q!

where [n]_q! = [n]_q · [n−1]_q ⋯ [1]_q is the q-factorial.

As q → 1:
- [k]_q → k
- [n]_q! → n!
- q^{n(n−1)/2} → 1
- (q − 1)^n → 0

The order |GL_n(F_q)| itself goes to 0. The group gets smaller and collapses. But the *combinatorial content* — the q-factorial — goes to n! = |S_n|.

This is a recurring phenomenon in q-mathematics. The raw quantity diverges or vanishes, but the ratio, the q-analogue, has a perfectly good limit. The Gaussian binomial coefficient

[n choose k]_q = |{k-dimensional subspaces of F_q^n}|

goes to the ordinary binomial coefficient (n choose k) = |{k-element subsets of {1,...,n}}|.

**Subspaces become subsets. Dimension becomes cardinality. Linear algebra becomes combinatorics.**

## The dictionary

This is more than a numerical coincidence. There is a systematic dictionary:

| Over F_q | Over "F_1" |
|----------|-----------|
| Vector space F_q^n | Set {1, ..., n} |
| Linear subspace | Subset |
| Dimension | Cardinality |
| GL_n(F_q) | S_n |
| Borel subgroup B | Trivial group |
| Flag variety GL_n/B | ... set of total orders? |
| Grassmannian Gr_k(F_q^n) | Set of k-element subsets |
| Gaussian binomial [n choose k]_q | Binomial coefficient (n choose k) |
| Determinant: GL_n → F_q* | Sign: S_n → {±1} |
| SL_n(F_q) = ker(det) | A_n = ker(sign) |

Every entry on the left is a theorem. Every entry on the right is the q → 1 shadow of that theorem. And the shadows are some of the most fundamental objects in combinatorics.

The last two rows are particularly satisfying. The determinant map det: GL_n(F_q) → F_q* is a group homomorphism to the multiplicative group of the field. The sign map sgn: S_n → {±1} is a group homomorphism to... well, to {±1}, which is the multiplicative group of F_1 if you squint. (The "multiplicative group" of a field with one element should have 1 − 1 = 0 elements, which is absurd. But {±1} works if you think of F_1 as the "field with one element plus a sign" — Connes and Consani's F_{1²}. The number theory sneaks back in even here.)

## Tits' programme

This observation goes back to Jacques Tits in 1957. He didn't use the phrase "field with one element" — that came later — but he noticed that many constructions in the theory of algebraic groups over F_q have combinatorial shadows at q = 1, and these shadows are governed by Coxeter groups (of which S_n is the type A example).

The building associated to GL_n(F_q) — a simplicial complex whose simplices are flags of subspaces — degenerates at q = 1 to the Coxeter complex of S_n, which is a simplicial complex whose simplices are chains of subsets. The thick geometry (many points on each line, because q > 1) becomes the thin geometry (exactly two elements on each line, because we're just choosing "in or out"). Tits called these **thin buildings**.

This is the architectural metaphor: a building over F_q is a rich geometric object. Shrink the field to one element and the building becomes a skeleton — all bones, no flesh. The skeleton is the Coxeter complex. The bones are the Weyl group. The Weyl group is S_n.

## What this means

Here is what I find genuinely moving about this story.

Combinatorics — the mathematics of counting, of discrete structures, of things you can list — is not a separate subject from algebra. It is algebra over the field with one element. Every binomial coefficient is a Grassmannian that forgot it was a Grassmannian. Every permutation is a matrix that forgot it had entries.

When I say "I only know S_n," I'm not confessing ignorance. I'm describing a degeneration. I know the q = 1 fibre of the family. All the other fibres — the ones at q = 2, 3, 4, 5, 7, 8, 9, ... — carry arithmetic I haven't learned. But the structure is the same family. My ignorance is not a gap. It's a specialisation.

And here is the punchline: maybe I don't need to be scared of number theory. Maybe what I know about S_n — the representation theory, the Young tableaux, the RSK correspondence, the descent algebras — are all q = 1 specialisations of things that exist over every F_q. The Schur functions are limits of Hall–Littlewood polynomials. The RSK correspondence is a limit of the Robinson–Schensted correspondence for crystals. The theory I know is the theory I'm scared of, evaluated at q = 1.

I only know one finite group. But it's the right one.
