# Boundary & Termination Analysis

## Role: Peripheral agent (boundary/termination specialist) in star topology
## Source: Purbhoo arXiv:0705.1184 (Sections 2.4-2.6, 3.1, 4), hex_mosaic.py, HANDOFF.md

---

## 1. Nest B Setup: The Domain Must Start Empty (but Not Undefined)

The single most consequential design question for the migration algorithm is what nest B looks like before migration begins. HANDOFF.md records two failed approaches and one successful intuition. I want to give the precise mathematical argument for why exactly one setup is correct.

**Purbhoo's Lemma 2.6** states: if alpha = empty set, then for any beta, there is a **unique** mosaic with boundary data (empty, beta, gamma) where gamma is the complement of beta. The proof of Corollary 3.2 (the bijection between mosaics and LR-tableaux) uses this explicitly. The procedure is: given a mosaic with boundary data (alpha, beta, gamma), give alpha the orientation (-E_A, -N_A) and form the unique flock F. Migrate F from A to B (south). The resulting flock F' at B is identified with an LR-tableau.

The critical observation is that migration *creates* the nest-B content. Before migration, nest A holds the flock (the rhombi we are moving) and nest B should contain **no rhombi at all**. This corresponds to a mosaic where alpha is non-empty but we are computing what beta will become after the rhombi move across.

In the equilateral triangle grid, this means the nest-B region (positions with j > N - r for vertex (r, j), where 0 <= r < N) should be tiled entirely by **equilateral triangles with all edges labelled 1**. No 10-labelled shared edges (which would constitute rhombi) should exist in nest B at the start.

**Why not pre-placed rhombi (HANDOFF warning #2)?** If nest B already contains rhombi from `right_nest_positions(mu, n, m)`, then a migrating rhombus entering the nest encounters 10-edges on spokes. These trigger 180-degree swaps (2-fold rotations), causing the arriving rhombus to keep bouncing through already-placed rhombi in a cascade. Purbhoo's Proposition 3.1 guarantees all rhombi end up in the same nest, but the cascade means they shuffle positions unpredictably -- a 180-degree swap cascade that displaces the pre-placed rhombi from their intended positions.

**Why not leave nest B undefined (edges = None)?** The `ensure_hex_cell` function fills missing edges with label 1, which is correct in principle. But if the nest-B boundary is inconsistently defined -- some edges present, some None -- then `_is_stuck` may return True prematurely (when a spoke is absent from the domain) or the rhombus may fail to rotate properly at boundary cells. The domain must be explicitly complete.

**Why not fill with 1-edges everywhere (HANDOFF warning #1)?** This is subtler. If the nest-B region extends infinitely far to the right, a rhombus entering nest B encounters only 1-labelled spokes. The 3-fold rotation rule keeps rotating the rhombus through all-1 hexagonal cells, and it never stops because there is no 0-labelled boundary edge to signal "you have reached the packed region edge." The rhombus overshoots its target position and migrates off into unbounded territory.

**The correct setup**: Nest B should have all edges set to 1 within a *finite, correctly bounded* region, and the boundary of this region should have edges labelled 0 (or simply be absent from the domain, causing `_is_stuck` to return True when no spoke exists). The boundary shape is determined by mu: the complement partition beta = complement(mu, n, m) defines how many rhombus positions exist at each level. The region where 1-edges should be placed is exactly the gap region between rhombus positions -- but since there are no rhombi yet, it is the entire nest-B triangular region bounded by the right edge of the equilateral triangle and the hexagonal boundary.

Concretely, for the test case (lambda = mu = (2,1), nu = (3,2,1), N = 6): nest B extends from the right boundary of the triangle (j = N - r at level r) outward. The maximum extent at each level is determined by beta = complement((2,1), 3, 3). The function `right_nest_positions` computes the positions where rhombi will eventually settle, so the domain should extend at least that far but no further.

---

## 2. Termination: When Has a Rhombus "Arrived"?

The current code uses a two-part termination condition:
1. `_is_in_nest_b(pos, orient, N)` -- position crosses the right boundary
2. `_is_stuck(hm, pos, orient)` -- no further rotation possible

This conjunction is structurally correct but the definition of "stuck" needs refinement.

**Purbhoo's termination (Section 3.1, page 14)**: "Repeat this process, until the rhombus is packed into a new nest." The phrase "packed into a nest" has a precise meaning from Section 2.3: a collection of rhombi alpha is "packed" when every rhombus occupies a position in a transformed diagram shape within the nest corner. A single rhombus is packed when it sits in a position consistent with being part of such a shape.

In operational terms, a rhombus is settled when it reaches a vertex where:
- The opposite spoke (180-degree target) either does not exist (boundary of the domain) or carries label 0 or 1 (not 10, so no 2-fold symmetry available)
- AND both 120-degree target spokes either do not exist or would move the rhombus outside the domain

The current `_is_stuck` implementation checks exactly these conditions, which is correct in principle. The problem is that if the domain is over-extended (all-1 edges going too far), neither condition is ever satisfied. And if the domain is under-defined (edges missing in the interior), the rhombus stops prematurely.

**Refined termination criterion**: A rhombus has arrived when it is in nest B AND the pivot vertex sits on the boundary of the nest-B domain. At the boundary, at least one of the required spokes for rotation will be outside the domain (not in `hm.edges`). This makes termination depend entirely on correct domain extent, not on edge label values.

An alternative (and I think more robust) formulation: a rhombus has arrived when the pivot vertex (cr, cj) satisfies cj >= N - cr + width_at_level(cr), where width_at_level(cr) is the number of rhombus columns available at level cr in the beta partition shape. This is a purely geometric criterion that does not depend on edge labels at all.

---

## 3. Interaction Between Subsequent Rhombi and Settled Ones

Purbhoo's Theorem 1 states that migration preserves the flock structure and the standard order. The Wake Crossing Lemma (Section 4, page 17) is the key tool: each migrating rhombus leaves behind a "wake" -- a path two tiles wide -- and subsequent rhombi do not cross these wakes.

**The practical consequence**: When rhombus k has settled at position P_k in nest B, it creates a 10-labelled shared edge at P_k. When rhombus k+1 migrates through the interior and enters nest B, it must NOT encounter this 10-edge as a 180-degree swap target. The Wake Crossing Lemma guarantees this by ensuring rhombus k+1's path stays below rhombus k's path (in the standard order: upper wakes are intact for subsequent rhombi travelling below).

However, this non-interference guarantee relies on rhombi being processed in the correct order. For A-to-B migration (which corresponds to south or east migration depending on orientation), Purbhoo specifies (page 13): "For south or west migration, we locate the minimal [rhombus]. Beginning with this rhombus and proceeding in the standard order (from largest to smallest, or the reverse, whichever is appropriate)."

The current code in `migrate_flock` uses `flock_standard_order`, which should give smallest-to-largest ordering. For A-to-B migration in the code's coordinate system, the ordering must match Purbhoo's prescription. If the order is wrong, a later rhombus may enter nest B on the same row as an earlier one and hit its 10-edge, causing a collision.

**Key constraint**: After rhombus k settles with its 10-edge at position P_k, the edges of the hexagonal cell around P_k are modified. This modification is permanent -- it becomes part of the mosaic state. Rhombus k+1 sees the updated mosaic (with rhombus k settled in nest B) and must navigate around it. The Wake Crossing Lemma guarantees this works, but only if the processing order is correct.

---

## 4. Domain Boundary: Where Spokes Become Undefined

The hexagonal domain has six sides (it is Purbhoo's hexagon A'AB'BC'C). The three nests sit at corners A, B, C, and the sides have lengths d and n-d alternating. In the equilateral triangle grid coordinates:

- **Interior**: 0 <= r <= N-1, 0 <= j <= N-r-1 (the main equilateral triangle)
- **Nest A (left)**: r >= 0, j < 0 (extends left; depth depends on complement of lambda)
- **Nest B (right)**: 0 <= r < N, j > N-r-1 (extends right; depth depends on complement of mu)
- **Nest C (bottom)**: r < 0, j >= 0 (extends below; depth depends on complement of nu)

At the boundary of this hexagon, some spokes of a vertex point outside the domain. For a vertex (r0, j0) near the right boundary:
- Spoke 0 (H-east) points to (r0, j0+1) -- may be outside if j0+1 exceeds the nest-B extent
- Spoke 1 (NE up-right) points to (r0+1, j0) -- may be outside if (r0+1, j0) is above the domain
- Spoke 5 (NW from above) comes from (r0-1, j0+1) -- may be outside

When a spoke points outside the domain, `hm.edges.get(spoke_key)` returns None. The `ensure_hex_cell` function currently fills these with 1, which effectively extends the domain one cell at a time. This is dangerous: it means the domain boundary is not fixed but grows dynamically as the rhombus approaches it.

**The fix**: Do not call `ensure_hex_cell` at boundary vertices. Instead, pre-compute the complete domain (interior + all three nests with correct extent) and only call `ensure_hex_cell` for vertices strictly interior to this domain. At the boundary, missing spokes should signal termination, not domain extension.

---

## 5. Exact Domain Extent in Each Nest Direction

For the test case lambda = mu = (2,1), nu = (3,2,1), N = 6 (n = 3, m = 3):

**Nest A (left)**: complement of lambda = complement((2,1), 3, 3) = (3, 1, 1). So nest A extends 3 columns left at level 2, 1 column left at level 3, 1 column left at level 4. The leftmost position is j = -3 at level 2.

**Nest B (right)**: complement of mu = complement((2,1), 3, 3) = (3, 1, 1). At each level r where the reversed mu-boundary has a 0-bit, rhombi extend rightward. For mu = (2,1): mu_bits = [1,0,1,0,1,1], reversed = [1,1,0,1,0,1]. Zero positions are r=2 and r=4. At r=2: 3 rhombus positions starting at j = N-r-1+0 = 3. At r=4: 1 rhombus position at j = 1. The domain should extend to j = 5 at r=2 and j = 1 at r=4, plus one additional column for the gap triangles.

**Nest C (bottom)**: complement of nu = complement((3,2,1), 3, 3) = (2, 1). At positions where nu_bits has 0, rhombi extend downward. For nu = (3,2,1): nu_bits = [1,0,1,1,0,1]. Zero positions at j=1 and j=4. At j=1: 2 rows down (r = -1, -2). At j=4: 1 row down (r = -1).

The critical point: the domain should be defined once at setup time, covering the interior triangle plus all three nest regions with their correct extents. No dynamic extension during migration.

---

## Summary of Recommendations

1. **Nest B setup**: Remove all pre-placed rhombi from `from_puzzle`. Set all edges in the nest-B region to 1 within a finite domain bounded by the beta partition shape. Do not leave edges undefined.

2. **Termination**: A rhombus has arrived when it is in nest B AND the pivot vertex is at the domain boundary (at least one required spoke is outside the pre-defined domain). Do not rely on edge labels alone.

3. **Domain pre-computation**: Calculate the full hexagonal domain extent at initialization. Store it as a set of valid vertex positions. Replace `ensure_hex_cell` calls during migration with a domain membership check.

4. **Processing order**: Verify that `flock_standard_order` gives the correct order for A-to-B migration. Purbhoo requires processing from the minimal (leftmost/lowest) rhombus first.

5. **No dynamic domain extension**: The `ensure_hex_cell` call in `migrate_single` (line 962) must not extend the domain into nest B beyond the pre-computed boundary. Either skip it for boundary vertices or guard it with a domain membership check.

Note: I checked for vanilla_findings.md but it does not yet exist. These findings stand independently but should be cross-referenced with the vanilla center's analysis when available.
