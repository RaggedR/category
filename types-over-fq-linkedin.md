What Happens to Types Over a Finite Field?

In Haskell, types form a semiring. You can add (sum types) and multiply (product types), but you can't subtract. Bool has 2 inhabitants. Maybe Bool has 3. There's no type with -1 inhabitants.

This is a problem if you want a Hopf algebra of types, because a Hopf algebra needs an antipode — a map S that "cancels" elements against their own decomposition. Cancellation requires negatives. A semiring doesn't have them.

But what if you change the ground field?

Over a finite field F_q, every element has an additive inverse. The semiring obstruction vanishes. And something remarkable happens: "types" become representations of quivers over F_q, and counting how they fit together gives you a genuine Hopf algebra with a non-trivial antipode.

A quiver is a directed graph. A representation assigns an F_q-vector space to each vertex and a linear map to each arrow. The Hall algebra of these representations has:

A product — given representations M and N, count the number of short exact sequences 0 → N → L → M → 0. The number of ways N can sit inside L with quotient M is the structure constant. This count is finite because we're over F_q.

A coproduct — dual to the product, decomposing a representation into sub and quotient.

A grading — by dimension vector. Connected at grade zero.

An antipode — forced by the connected graded theorem, the same mechanism that gives Connes-Kreimer its antipode on rooted trees.

This is a real Hopf algebra built from "types over F_q." The types are representations. The algebra structure comes from how they extend and decompose. The antipode computes inclusion-exclusion over the lattice of subrepresentations.

Ringel proved that this Hall algebra is isomorphic to the positive part of the corresponding quantum group U_q(g). The q in "quantum group" is the size of the finite field. It was there all along — not a deformation parameter, but a field size.

So:

Types over the natural numbers — a semiring. No negatives, no antipode, no Hopf algebra.

Types over F_q — a Hall algebra. Genuine Hopf algebra. Isomorphic to half a quantum group.

Take q → 1 — recover the classical universal enveloping algebra of a Lie algebra.

The algebra of types doesn't admit an antipode over N, where Haskell lives. But move to F_q, and the antipode isn't just possible — it's a quantum group, one of the deepest structures in modern mathematics.

The semiring obstruction was real. It was also an artefact of the ground field. Change the field, and types don't just gain negatives. They gain an antipode, a quantum group, and a connection to representation theory that nobody would have guessed by staring at Haskell's type system.

Sometimes the interesting structure isn't in the object. It's in what happens when you change where the object lives.
