# Coordinate Translation: Purbhoo's Mosaic Geometry to hex_mosaic.py's Equilateral Grid

## 1. The Two Geometries

Purbhoo's mosaics (Section 2.4) tile a hexagon A'AB'BC'C with three tile types:
- (a) equilateral triangle, side 1
- (b) square, side 1
- (c) rhombus, side 1, angles 30 and 150 degrees

The hexagon has 150-degree angles at the nest corners A, B, C and 90-degree angles at the intermediate vertices A', B', C'. The sides alternate between length d (nest sides: A'A, B'B, C'C) and n-d (interior sides: AB', BC', CA').

hex_mosaic.py uses a different geometry: a standard equilateral triangle grid where every face is an equilateral triangle (either "up" or "down"), every edge connects two vertices, and every vertex has exactly 6 spokes. The three edge orientations are H (horizontal), NE (northeast), and NW (northwest). There are no squares and no 30/150-degree rhombi in this grid -- only equilateral triangles and their pairings into 60/120-degree rhombi.

The fundamental question is: how does the mosaic's geometry, with its squares and thin rhombi, map onto this equilateral grid where the only tile shapes are equilateral triangles?

## 2. The Bijection Between Mosaics and Puzzles

Purbhoo describes the bijection physically (page 11): take a mosaic, remove all rhombi, and imagine the squares have hinges at their corners. Grab A', B', C' and pull outward until the boundary straightens into an equilateral triangle of side n. The squares flex into 60/120-degree rhombi. Some line segments rotate clockwise (labeled 0), others counterclockwise (labeled 1). This produces the puzzle.

This bijection tells us how tiles transform:
- **Mosaic triangle** (a) -> remains a single triangle face in the grid, with all edges carrying the same label (all 0 or all 1).
- **Mosaic square** (b) -> flexes into a 60/120-degree rhombus in the grid, composed of one up-face and one down-face sharing an internal edge. The shared edge carries label 10 in the edge-label scheme. From the DELTA/NABLA analysis, this is the H-oriented rhombus: the shared edge is horizontal, edge labels are (10, 1, 0) on the up-face and (10, 0, 1) on the down-face.
- **Mosaic rhombus** (c) -> becomes a 60/120-degree rhombus in the grid, but with a different internal edge orientation. Since 30/150-degree rhombi come in three orientations (one per nest), each nest's rhombi use a distinct shared-edge orientation.

## 3. Tile Type to Grid Edge Label Mapping

In hex_mosaic.py, a rhombus is identified by which edge carries the label 10 (the "shared edge" between its two triangular faces). The three possibilities:

| Grid edge type | Shared edge formula | Tile identity in Purbhoo |
|---|---|---|
| NW (label 10 on NW edge) | ((r, j+1), (r+1, j)) | Nest A (left) rhombus -- a 30/150 mosaic rhombus (c) |
| H (label 10 on H edge) | ((r+1, j), (r+1, j+1)) | Mosaic square (b) when in the interior; nest B (right) rhombus when in the right nest |
| NE (label 10 on NE edge) | ((r, j+1), (r+1, j+1)) | Nest C (bottom) rhombus -- a 30/150 mosaic rhombus (c) |

This is confirmed by the code's nest placement functions:
- `left_nest_positions` (lines 45-65): places NW-oriented rhombi at j < 0, with `_nw_edge(r, j) = edge_key((r, j+1), (r+1, j))`.
- `right_nest_positions` (lines 68-93): places H-oriented rhombi at j > N-r, with `_h_edge(r, j) = edge_key((r+1, j), (r+1, j+1))`.
- `bottom_nest_positions` (lines 96-119): places NE-oriented rhombi at r < 0, with `_ne_edge(r, j) = edge_key((r, j+1), (r+1, j+1))`.

The critical insight is that Purbhoo's three rhombus orientations in the mosaic (one per nest) correspond to the three edge orientations in the equilateral grid. Each nest uses exactly one orientation, and these are NW, H, NE for nests A, B, C respectively. The mosaic squares (which are interior tiles, not nest tiles) also map to H-oriented rhombi in the grid -- they share the H orientation with nest B, which is natural because nest B is the "right" nest and squares flex horizontally when the boundary is straightened.

## 4. Edge Orientations and Purbhoo's "Horizontal" and "Vertical"

Purbhoo's direction rule (page 13) says: for case (b), rotate 120 degrees "in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

The critical warning on page 4: "all directions in this paper are used in a relative, rather than in an absolute sense." For A-to-B migration, the diagrams should be rotated so the migration direction points right.

In the hexagon A'AB'BC'C (Figure 5), the orientation vectors at each nest are:
- Nest A: (E_A, N_A) parallel to sides AA' and AB' respectively
- Nest B: (E_B, N_B) parallel to sides BB' and BC' respectively
- Nest C: (E_C, N_C) parallel to sides CC' and CA' respectively

For A-to-B migration (south in Purbhoo's convention, since the flock at A migrates toward B), we orient the page so the migration direction points right. "Horizontal" and "vertical" are defined relative to this oriented view.

Now, in hex_mosaic.py's equilateral grid, the three edge orientations have geometric directions:
- **H edges**: truly horizontal, connecting (r, j) to (r, j+1) -- these run left-right
- **NE edges**: at 60 degrees from horizontal, connecting (r, j) to (r+1, j) -- these run from lower-left to upper-right
- **NW edges**: at 120 degrees from horizontal, connecting (r, j+1) to (r+1, j) -- these run from lower-right to upper-left

For A-to-B migration in the grid, the rhombus moves from left (j < 0) toward right (j > N-r). The natural "rightward" direction in the grid is increasing j, which is the H direction. So:
- **H = "horizontal"** in Purbhoo's sense (parallel to the migration direction)
- **NE = "vertical"** in Purbhoo's sense (perpendicular to migration, going "up" in the oriented view -- 60 degrees from horizontal in the grid corresponds to the direction perpendicular to the rightward migration axis)
- **NW = neither horizontal nor vertical** -- this is the orientation that should be ruled out by the direction rule

This is exactly the interpretation that the current code uses (line 976-979):
```python
# For A→B migration: H = horizontal, NE = vertical, NW = forbidden.
# NE spokes (1, 4) → CW; all others (H: 0,3 and NW: 2,5) → CCW.
direction = 1 if my_spoke in (1, 4) else -1
```

## 5. The Ambiguity at NW Spokes

The direction rule table from CONTEXT.md reveals the problem. When the current spoke is NW (spokes 2 or 5), both rotation directions produce a "good" orientation:

| my_spoke | orient | CW result | CCW result |
|---|---|---|---|
| 2 | NW | spoke 4 (NE) | spoke 0 (H) |
| 5 | NW | spoke 1 (NE) | spoke 3 (H) |

Both NE and H are "allowed" by the rule (both are horizontal-or-vertical). The current code resolves this by choosing CCW for NW spokes (giving H), but it is not obvious from the paper alone that this is correct.

However, there is a structural argument. The rhombus at a NW spoke is currently oriented as a nest-A rhombus. During A-to-B migration, it needs to transition toward nest-B orientation (H). The most natural step in the orientation cycle is NW -> H (one 120-degree rotation toward the target nest), not NW -> NE (which overshoots past H toward nest-C orientation). So CCW (giving H) is the geometrically correct choice for NW spokes.

This gives the complete direction rule for A-to-B migration:
- **spoke 0 (H)**: CCW -> spoke 4 (NE). Correct: H is already "horizontal," rotation advances toward "vertical."
- **spoke 1 (NE)**: CW -> spoke 3 (H). Correct: NE is "vertical," rotation brings it back to "horizontal."
- **spoke 2 (NW)**: CCW -> spoke 0 (H). Correct: NW advances to H (toward target nest).
- **spoke 3 (H)**: CCW -> spoke 1 (NE). Same as spoke 0 by symmetry.
- **spoke 4 (NE)**: CW -> spoke 0 (H). Same as spoke 1 by symmetry.
- **spoke 5 (NW)**: CCW -> spoke 3 (H). Same as spoke 2 by symmetry.

The orientation cycle is: NW -> H -> NE -> NW -> ..., always advancing in the CCW direction around the spoke wheel. But the *rotation* direction alternates: CCW for NW and H spokes, CW for NE spokes.

## 6. Nest Extent in Grid Coordinates

Each nest is a cone of angle 150 degrees at its corner vertex.

**Nest A (left)**: Corner at the left vertex of the triangle. Rhombi have NW orientation. In grid coordinates, nest A occupies positions (r, j) with j < 0 and 0 <= r < N. The extent is governed by the partition alpha = complement(lambda): row r has alpha_i rhombi extending leftward from j = -1 to j = -alpha_i, where the row indices are determined by the 0-bits in the lambda boundary string. Maximum extent: j = -max(alpha_i) = -(n-d) when alpha has a full row.

**Nest B (right)**: Corner at the right vertex. Rhombi have H orientation. In grid coordinates, nest B occupies positions (r, j) with j >= N-r-1 and 0 <= r < N. The right boundary at level r is at column N-r, and rhombi extend rightward from j = N-r-1 to j = N-r-1+beta_i-1, where beta = complement(mu) and row indices come from the reversed mu boundary bits. Maximum extent: j = N + (n-d) - 1.

**Nest C (bottom)**: Corner at the bottom vertex. Rhombi have NE orientation. In grid coordinates, nest C occupies positions (r, j) with r < 0 and 0 <= j < N. Rhombi extend downward from r = -1 to r = -gamma_i, where gamma = complement(nu) and column indices come from the nu boundary bits. Maximum extent: r = -(n-d).

## 7. Summary and Remaining Questions

The coordinate translation is now well-defined:

1. Purbhoo's three mosaic tile types map to equilateral grid structures as: triangle -> single face, square -> H-oriented rhombus (label 10 on horizontal edge), 30/150-rhombus -> NW or NE-oriented rhombus depending on nest.

2. The grid's three edge orientations H, NE, NW correspond to Purbhoo's relative "horizontal," "vertical," and "forbidden" (the one ruled out by the direction rule) for A-to-B migration.

3. Each nest uses one rhombus orientation: A=NW, B=H, C=NE.

4. Nest extents are determined by complement partitions within cones of known grid-coordinate boundaries.

The direction rule as currently implemented (CCW for H and NW spokes, CW for NE spokes) appears correct by the structural argument that the orientation should advance NW -> H -> NE in the same cyclic direction, always preferring the rotation that moves toward (not past) the target orientation. If collisions persist, the bug is likely elsewhere -- perhaps in the pivot vertex selection, the spoke-to-position mapping, or the order in which rhombi are processed.

*No vanilla_findings.md was available at the time of writing.*
