My Favourite Concept from Category Theory

Hopf algebras.

A Hopf algebra is five maps on one object, constrained by commutative diagrams until only one consistent structure remains. Let me walk through it.

Start with a vector space H over a field k.

An algebra gives you a way to combine: multiplication m: H ⊗ H → H and unit η: k → H, satisfying associativity:

m(m ⊗ id) = m(id ⊗ m).

A coalgebra gives you a way to decompose: comultiplication Δ: H → H ⊗ H and counit ε: H → k, satisfying the same axiom with arrows reversed:

(Δ ⊗ id)Δ = (id ⊗ Δ)Δ.

The maps: Δ and ε must be algebra morphisms, or equivalently, m and η must be coalgebra morphisms. This compatibility gives you a bialgebra.

A Hopf algebra adds one more map: the antipode S: H → H, satisfying:

m(S ⊗ id)Δ = η ∘ ε = m(id ⊗ S)Δ.

The same axiom, instantiated in different Hopf algebras, produces remarkably different-looking maps. In every graded Hopf algebra the map S is Mobius inversion on the poset determined by the coproduct:

Polynomials k[x] — the poset is trivial. S(x) = -x. Mobius inversion reduces to sign alternation.

Group algebras k[G] — not graded, the map S is group inversion. S(g) = g⁻¹.

Symmetric functions Sym — the poset is partitions ordered by refinement. S(s_λ) = (-1)^|λ| s_{λ'}, where λ' is the conjugate partition. Mobius inversion transposes Young diagrams.

Other interesting Hopf algebras:
- QSym (commutative) / NSym (non-commutative) — dual pair
- Connes-Kreimer H_R (commutative, basis = rooted forests) / Grossman-Larson H_GL (non-commutative, basis = rooted trees (grading shift by one)) — dual pair
- Foissy H_PR (non-commutative, basis = planar rooted forests, dimensions = Catalan numbers) — self-dual

These Hopf algebras are related to each other in beautiful ways that I don't have space to go into here.

Now for my "open question".

Haskell's type system already borrows heavily from category theory. Types form a category. Functors, monads, products, and coproducts are all first-class concepts. Product types are categorical products. Sum types are coproducts. The category of types is cartesian closed.

A Hopf algebra needs more: a vector space structure, a tensor product, a compatible algebra and coalgebra structure, and an antipode. The product and coproduct of types give you something — but are they compatible in the bialgebra sense?

What would it mean for an algebra of types to admit an antipode?

I would love to know.
