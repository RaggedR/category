Every cell in the Game of Life is solving the same problem: "what's happening around me?"

The rules are simple. Each cell on an infinite grid is alive or dead. Every generation, it looks at its 8 neighbors: survive with 2-3 neighbors, be born with exactly 3, die otherwise. From this, you get gliders, oscillators, universal computation.

But notice: the cell doesn't just receive a value. It receives a perspective — the whole grid, seen from where it sits. It peeks in every direction, reads its neighborhood, and collapses that view into a single decision.

Every cell does this simultaneously. Same rule, different vantage point.

Mathematicians call this a comonad. The life rule has the shape Grid a → a: consume an ambient context, produce a local result. This shape is called a co-Kleisli arrow.

Now the real question: what does a cell look like after two generations?

You can't get this by peeking at your neighborhood once. To know your state at generation 2, you need your neighbors' states at generation 1. But their gen-1 states each depended on their neighbors at generation 0. You need information aggregated from two hops away.

Co-Kleisli composition solves this in two steps:

extend(life_rule): at every position, set up a perspective (the whole grid focused there), then apply life_rule from that perspective. Every cell computes its own next state simultaneously. This builds generation 1.

life_rule: apply the rule one more time at your position, now reading from generation 1.

The formula is: life_rule <=< life_rule = life_rule ∘ extend(life_rule)

And the result is still Grid a → a. Same shape. You can compose again. One life_rule peeks at a 3×3 block (9 cells). Two compositions reach a 5×5 block (25 cells). N compositions reach (2N+1)² cells. Each composition widens the effective neighborhood without changing the interface.

What surprised me is how close this is to a concept from algebraic geometry: sheaves. A sheaf assigns local data to each point of a space, with the requirement that overlapping regions agree — local observations, globally coherent. The "set up a perspective at every position" step inside extend does exactly this: it builds a "section" at every point (the grid focused there), and the comonad laws guarantee these sections are compatible.

The Game of Life is a sheaf of focused grids, with co-Kleisli composition layering new sheaves on top of old ones.

Three things that don't usually appear together: a toy from recreational computing, a design pattern from functional programming, a concept from algebraic geometry. Same idea, different clothes.
