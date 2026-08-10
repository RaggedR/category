# Agent 2: Coordinate Translator — Mapping Purbhoo's Mosaic Geometry to the Equilateral Grid

## 1. The Two Coordinate Systems

Purbhoo's mosaics live in a hexagonal domain A'AB'BC'C with angles 150 at A,B,C and 90 at A',B',C'. The hex_mosaic.py code works in an equilateral triangle grid with vertex coordinates (r, j), where r is the level (0 = bottom) and j is the horizontal position. These are fundamentally different geometric embeddings of the same combinatorial structure. The key to fixing the direction rule is understanding how they relate.

## 2. The Hexagon A'AB'BC'C and Its Embedding in the Grid

Purbhoo's hexagon (Section 2.4, Figure 5) has:
- Side lengths: A'A = B'B = C'C = d, and AB' = BC' = CA' = n-d
- Three nests at corners A, B, C (the 150-degree corners)
- Three primed vertices A', B', C' (the 90-degree corners)

In hex_mosaic.py, the hexagonal domain is embedded so that:
- The **interior** is the equilateral triangle with r=0..N-1, j=0..N-r-1 (same as the puzzle)
- **Nest A** (left) extends into j < 0, i.e. to the LEFT of the triangle
- **Nest B** (right) extends into j > N-r, i.e. to the RIGHT of the triangle
- **Nest C** (bottom) extends into r < 0, i.e. BELOW the triangle

The hexagon vertices map approximately as follows:
- **A** = top-left corner of the interior triangle ~ vertex (N, 0) (the apex, where left and right boundaries meet for nest A)
- **B** = top-right ~ far right boundary region
- **C** = bottom ~ base of the triangle
- **A'** = bottom-left boundary (where left nest meets bottom nest)
- **B'** = along the right boundary (where interior meets right nest)
- **C'** = along the bottom boundary (where interior meets bottom nest)

The crucial point: **the Purbhoo hexagon is not drawn as an equilateral triangle.** The straightening bijection (Section 2.4, page 11) describes physically pulling the mosaic by its primed vertices to flatten it into an equilateral triangle + three nest regions. The squares in the mosaic become rhombi in the equilateral grid, and vice versa.

## 3. Tile Type Correspondence

Purbhoo defines three tile types (Section 2.4, page 11):
- **(a) Equilateral triangle**, side 1
- **(b) Square**, side 1
- **(c) Rhombus**, side 1, angles 30 and 150 degrees

In the equilateral triangle grid of hex_mosaic.py, ALL tiles are composed of pairs of adjacent up/down triangular faces. The edge labels {0, 1, 10} determine the tile type. From the established findings in the memory file and the code:

| Purbhoo tile | Grid representation | Shared edge label | Edge triple (DELTA/NABLA) |
|---|---|---|---|
| Triangle (a) | Single triangular face | n/a (all edges 0 or all 1) | (0,0,0) or (1,1,1) |
| Square (b) | Two faces sharing an **H edge** with label 10 | 10 on H | (10,1,0) or (10,0,1) depending on up/down |
| Rhombus (c), 30/150 | Two faces sharing a **NW** or **NE edge** with label 10 | 10 on NW or NE | (1,0,10) or (0,10,1) etc. |

This is the critical mapping: **Purbhoo's squares correspond to H-oriented rhombi in the grid, and Purbhoo's 30/150 rhombi correspond to NW-oriented or NE-oriented rhombi in the grid.**

The reason: when the hexagonal mosaic is "straightened" into the equilateral triangle grid, the right angles of the squares get distorted into the 60/120 angles of H-rhombi, while the 30/150 rhombi in the mosaic retain their non-equilateral character as NW or NE pairs.

## 4. "Horizontal" and "Vertical" in Purbhoo vs. Grid Edge Orientations

This is the heart of the bug. Purbhoo says (page 13): "rotate 120 degrees, in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

Purbhoo explicitly warns (page 4): "all directions in this paper are used in a relative, rather than in an absolute sense." For migration from A to B, the diagram should be rotated so the migration direction points **to the right** (page 13).

In the hex_mosaic.py code, migration from nest A (left, j<0) to nest B (right, j>N-r) means moving rightward in the j direction. This aligns with Purbhoo's convention.

Now, what are "horizontal" and "vertical" relative to the A-to-B migration direction?

**In Purbhoo's mosaic geometry** (the physical hexagon, NOT the equilateral grid):
- **Horizontal** = parallel to the migration direction (east-west in the diagram). This is the direction along the **AA'** or **BB'** edges of the hexagon.
- **Vertical** = perpendicular to the migration direction (north-south in the diagram). This is the direction along the **AB'** interior boundary.

**In the equilateral grid:**
The three edge orientations are:
- **H** (horizontal): connects (r,j) to (r,j+1) -- truly horizontal in the grid
- **NE** (northeast /): connects (r,j) to (r+1,j) -- slopes upward-left at 60 degrees
- **NW** (northwest \\): connects (r,j+1) to (r+1,j) -- slopes upward-right at 60 degrees

The hexagon edge directions from Section 2.4 define orientation vectors:
- E_A, N_A are in the directions of AA', AB' respectively
- E_B, N_B are in the directions of BB', BC' respectively

For A-to-B migration, the migration direction is roughly from A toward B. Looking at the hexagon geometry (Figure 5), the side AB' has length n-d and runs between the two nests. The edges AA' and BB' (length d) are the nest-internal sides.

**The key correspondence is:**

| Purbhoo direction concept | Hexagon side | Grid edge orientation |
|---|---|---|
| "Horizontal" (parallel to migration) | Along AB' (interior, connecting nests) | **H** edges in the grid |
| "Vertical" (perpendicular to migration) | Along AA' or BB' (nest sides) | **NE** edges in the grid |
| Neither | The third direction | **NW** edges in the grid |

The reasoning: In the equilateral grid, H edges run horizontally (constant r), and migration from left (small j) to right (large j) is also horizontal. So H edges are "horizontal" in Purbhoo's sense. NE edges connect different levels (different r values), running perpendicular to the horizontal migration flow -- these are "vertical." NW edges are the third direction -- neither horizontal nor vertical.

## 5. The Direction Rule Translated

Purbhoo's rule: "rotate 120 degrees in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

Translated to the grid: **rotate 120 degrees in whichever way causes the rhombus's shared edge to become H or NE orientation (NOT NW).**

This matches the interpretation already in the code comment at line 711:
```
# In mosaic geometry: H = horizontal, NE = vertical, NW = neither.
# After rotation, shared edge must be H or NE, not NW.
```

And it produces the `_MIGRATION_SHIFT` table at line 714:
```python
_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}
```

Let me verify this against the spoke-to-orient table from CONTEXT.md:

| my_spoke | orient | CW -> spoke (orient) | CCW -> spoke (orient) | Rule: pick non-NW | Result |
|---|---|---|---|---|---|
| 0 | H | 2 (NW) | 4 (NE) | NE wins | **CCW**, shift=2 |
| 1 | NE | 3 (H) | 5 (NW) | H wins | **CW**, shift=4 |
| 2 | NW | 4 (NE) | 0 (H) | **BOTH non-NW** | AMBIGUOUS |
| 3 | H | 5 (NW) | 1 (NE) | NE wins | **CCW**, shift=2 |
| 4 | NE | 0 (H) | 2 (NW) | H wins | **CW**, shift=4 |
| 5 | NW | 1 (NE) | 3 (H) | **BOTH non-NW** | AMBIGUOUS |

**Spokes 2 and 5 (the NW-oriented ones) are ambiguous** -- both CW and CCW produce a non-NW result. The `_MIGRATION_SHIFT` table assigns shift=2 (CCW) to both, choosing H over NE. But this needs justification.

## 6. Resolving the NW Spoke Ambiguity

For spoke 2 (NW): CW gives spoke 4 (NE), CCW gives spoke 0 (H).
For spoke 5 (NW): CW gives spoke 1 (NE), CCW gives spoke 3 (H).

In both ambiguous cases, one option gives H and the other gives NE. Since Purbhoo says "horizontal or vertical" (not "vertical or horizontal"), and since these words might have a preference order, one possibility is:
- **Prefer H** (horizontal) when both are available, because horizontal is listed first.

But this is speculative. A more principled resolution comes from the mosaic geometry:

When a NW-oriented rhombus sits in the mosaic, it has edges going in the "third" direction (neither horizontal nor vertical). After a 120-degree rotation, it needs to become either horizontal or vertical. The question is which 120-degree rotation advances the rhombus in the migration direction (rightward/toward nest B).

**Geometric argument:** For A-to-B migration (rightward), the rhombus should generally advance its j-coordinate. CCW rotation (shift=2) on NW spokes produces H orientation, which keeps the rhombus at the same level and shifts it horizontally -- this is the "correct" direction because it advances the rhombus toward nest B along the migration path. CW rotation would produce NE orientation, which changes the level -- a "vertical" move that doesn't directly advance migration.

However, this geometric argument may fail in specific configurations. The position-based heuristic in `migrate_single` was trying to capture this but gets the tiebreaks wrong.

## 7. Nest Positions in Grid Coordinates

From the code:

**Nest A (left):** Rhombi at positions (r, j) with j < 0. These are NW-oriented rhombi: each is a pair ('up', r, j) + ('down', r, j) sharing the NW edge at ((r, j+1), (r+1, j)). The grid levels (r values) correspond to 0-bits in the lambda boundary string.

**Nest B (right):** Rhombi at positions (r, j) with j > N-r-1. These are H-oriented rhombi: each is ('down', r, j) + ('up', r+1, j) sharing the H edge at ((r+1, j), (r+1, j+1)). The levels come from 0-bits in the reversed mu boundary string.

**Nest C (bottom):** Rhombi at positions (r, j) with r < 0. These are NE-oriented rhombi: each is ('down', r, j) + ('up', r, j+1) sharing the NE edge at ((r, j+1), (r+1, j+1)).

This is significant: **each nest uses a different rhombus orientation!** Nest A packs NW, nest B packs H, nest C packs NE. This means that during migration from A to B, rhombi must transition from NW orientation to H orientation. The intermediate orientations they pass through during migration are exactly what the direction rule governs.

## 8. Concrete Translation Table

| Purbhoo concept | hex_mosaic.py equivalent |
|---|---|
| Hexagon A'AB'BC'C | Triangle grid (r,j) + three nest extensions |
| Nest A (150-degree corner) | Left region, j < 0, NW-oriented rhombi |
| Nest B (150-degree corner) | Right region, j > N-r, H-oriented rhombi |
| Nest C (150-degree corner) | Bottom region, r < 0, NE-oriented rhombi |
| Triangle tile (a) | Single face with labels (0,0,0) or (1,1,1) |
| Square tile (b) | H-oriented rhombus (10 on H edge) |
| 30/150 rhombus tile (c) | NW or NE-oriented rhombus (10 on NW or NE edge) |
| "Horizontal" edge | H edge in grid: ((r,j),(r,j+1)) |
| "Vertical" edge | NE edge in grid: ((r,j),(r+1,j)) |
| Neither horiz nor vert | NW edge in grid: ((r,j+1),(r+1,j)) |
| E_A (east at nest A) | Direction of decreasing j (into the nest) |
| N_A (north at nest A) | Direction along the left boundary |
| E_B (east at nest B) | Direction of increasing j (into the nest) |
| N_B (north at nest B) | Direction along the right boundary |
| Migration A to B | Increasing j, moving right through interior |
| Flock ordering (standard) | By (label asc, j desc, r desc) |
| 2-fold symmetry (case a) | Opposite spoke has label 10 -> rotate 180 |
| 3-fold symmetry (case b) | Opposite spoke != 10 -> rotate 120 |
| Direction rule | Result orientation should be H or NE, not NW |
| Pivot vertex | Leading vertex toward nest B (higher j, then higher r) |

## 9. Cross-Reference with Agent 1 (Paper Exegete)

Agent 1's findings (agent1_findings.md) converge strongly with this analysis on all major points:

1. **Agreement on direction mapping:** Agent 1 independently concludes H = "horizontal" (parallel to migration), NE = "vertical" (perpendicular), NW = "oblique" (the forbidden direction). This matches Section 4 above exactly.

2. **Agreement on the bug:** Agent 1 confirms that `migrate_single` uses `direction_score` (a position-based heuristic) instead of the `_MIGRATION_SHIFT` table that already exists in the code (line 714). The table encodes the correct orientation-based rule; it just isn't called.

3. **Agreement on NW ambiguity:** Agent 1 reaches the same conclusion that spokes 2 and 5 are the ambiguous cases. Both CW (giving NE/"vertical") and CCW (giving H/"horizontal") satisfy the "not NW" constraint.

4. **Agent 1's resolution of the ambiguity:** Agent 1 cites Purbhoo's remark (page 13-14) that "the rhombus moves roughly to the right (although sometimes it moves strictly vertically)." This implies the default for NW spokes is CCW (producing H, i.e., horizontal/rightward motion), with CW (producing NE, vertical motion) as the exception. This supports the `_MIGRATION_SHIFT` table's choice of shift=2 (CCW) for spokes 2 and 5.

5. **Remaining open question:** Agent 1 raises whether the NW-spoke choice might be context-dependent (always CCW, or sometimes CW near boundaries). The `_MIGRATION_SHIFT` static table assumes always CCW. If the table is correct, this should be verifiable by tracing Figure 2.

## 10. The Unresolved Question

Two concerns remain:

**First**, does "some edge of the rhombus" in Purbhoo's rule mean the SHARED edge (the 10-labeled spoke), or ANY of the four rhombus edges? Agent 1 analyzed this carefully and concluded (via a contradiction argument) that it must refer to the shared edge specifically. If it meant any edge, the rule would be trivially satisfied by both rotations (since a rhombus always occupies exactly two of three directions, and there is only one "oblique" direction, at least one post-rotation direction is always non-oblique). So "some edge" = the shared/spoke edge.

**Second**, is the `_MIGRATION_SHIFT` table correct as-is, or does the NW case need context-sensitivity? The simplest hypothesis is that the static table is correct and the bug is purely that `migrate_single` doesn't use it. This is testable: replace the `direction_score` block with a `_MIGRATION_SHIFT[my_spoke]` lookup and run the test case (lambda=mu=(2,1), nu=(3,2,1), N=6). If all 6 rhombi land at 6 distinct positions, the table is correct. If collisions persist, the NW case needs refinement.
