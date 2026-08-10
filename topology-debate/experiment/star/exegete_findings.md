# Exegete Findings: The Exact Direction Rule for 120-Degree Rotations in Purbhoo's Migration

## Source Material
Purbhoo, "Puzzles, tableaux, and mosaics" (arXiv:0705.1184), primarily Section 3.1 (pages 12-14), with supporting material from Sections 2.3-2.4 (pages 8-11) and Section 4 (pages 17-18). Figure 2 (page 5) provides the only worked example.

## Note on Vanilla Findings
The vanilla center agent's findings file does not yet exist. This analysis proceeds from the paper alone.

---

## 1. The Exact Statement of the Direction Rule

The critical passage appears on page 13, in the description of what happens when a rhombus encounters a hexagon with 3-fold symmetry during migration:

> "If we see (b), the hexagon has a 3-fold rotational symmetry, so we rotate the tiling 120 degrees, **in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.** (This restriction rules out one of the two possible 120 degree rotations.)"

Purbhoo claims this is a complete disambiguation: of the two possible 120-degree rotations (clockwise and counterclockwise), exactly one produces an edge that is "horizontal or vertical." The other does not. The rule is deterministic.

But the terms "horizontal" and "vertical" are treacherous here. Let me trace exactly what they mean.

---

## 2. What "Horizontal or Vertical" Means in Purbhoo's Coordinate System

### 2.1 The Relativity Warning

On page 4, Purbhoo issues an explicit warning:

> "We would also like to warn the reader that all directions in this paper are used in a relative, rather than in an absolute sense. When we speak of 'up', 'down', 'left' and 'right', these are defined relative to a chosen 'right'; hence the reader may sometimes need to rotate the page."

Furthermore: "Compass directions 'north', 'east', 'south' and 'west', are used even more loosely. These are defined relative to a chosen basis for R^2, which may not be orthogonal or even have the expected orientation."

This means "horizontal" and "vertical" in the direction rule are NOT absolute page directions. They are defined relative to the migration direction.

### 2.2 The Migration Convention

On page 13, immediately before the direction rule:

> "In describing migration we will assume that our diagrams are turned so that the direction of migration is pointing directly to our right."

So the convention is: **rotate the page so migration points right**. Under this convention, "horizontal" means parallel to the migration direction, and "vertical" means perpendicular to it.

### 2.3 The Hexagonal Coordinate System

From Section 2.4 (page 10): The hexagon A'AB'BC'C has nests at corners A, B, C. The orientation vectors are defined by the hexagon edges:

- E_A, N_A are directions of edges AA', AB' respectively
- E_B, N_B are directions of edges BB', BC' respectively
- E_C, N_C are directions of edges CC', CA' respectively

These are NOT orthogonal pairs. They meet at 150 degrees (the angles at A, B, C). The angles at A', B', C' are 90 degrees.

### 2.4 Resolving the Semantics

For migration from A to B, the compass direction is one of four possibilities: north = N, east = E, south = -N, west = -E. The paper specifies (page 13) that the direction must be "consistent with the orientation of the flock." For nest A with orientation (E_A, N_A), the available migration directions are north (N_A) or east (E_A) if the nest opens to the northeast, and south (-N_A) or west (-E_A) otherwise.

Now, "horizontal" means parallel to the migration direction (right on the rotated page), and "vertical" means perpendicular to the migration direction (up/down on the rotated page).

**Key insight**: In the mosaic's equilateral triangle grid, there are three edge orientations. Once we fix a migration direction and orient it rightward, the three edge families become:
- One family roughly parallel to the migration direction (this is "horizontal")
- One family roughly perpendicular to the migration direction (this is "vertical")
- One family at an oblique angle to the migration direction (this is NEITHER horizontal nor vertical)

The direction rule says: choose the 120-degree rotation that makes the rhombus's shared edge become horizontal or vertical. The rotation that would make it oblique is the wrong choice.

---

## 3. Mapping to the hex_mosaic.py Coordinate System

### 3.1 The Three Edge Orientations

In the equilateral triangle grid used by `hex_mosaic.py`, edges come in three orientations:
- **H** (horizontal): edges connecting vertices at the same row
- **NE** (northeast): edges going up-right
- **NW** (northwest): edges going up-left

### 3.2 Migration Direction and the "Oblique" Orientation

For A-to-B migration, the migration direction depends on which edges of the hexagon connect A to B. The CONTEXT.md states the migration is coded as direction 'nw' in the implementation, meaning the rhombi move from the left nest toward the right nest.

The crucial question is: which of {H, NE, NW} is the "oblique" orientation -- the one that Purbhoo's rule EXCLUDES?

Consider the geometry. When the page is rotated so migration points right:
- If migration runs along the E_A direction, then edges parallel to E_A are "horizontal"
- Edges perpendicular to E_A (in the 150-degree sense) are "vertical"
- The third edge family, at 60 degrees to both, is "oblique"

For the standard hexagon placement (A at top-left, B at bottom, C at right -- see Figure 5), migration from A to B goes roughly south. Rotating the page so this points right means a 90-degree clockwise rotation of the diagram. Under this rotation:
- The H edges become roughly vertical
- The NE edges become roughly horizontal (parallel to migration)
- The NW edges become oblique (at 60 degrees to migration)

**Therefore, for A-to-B migration, NW is the excluded orientation.** The direction rule says: choose the rotation that gives H or NE, not NW.

### 3.3 The Ambiguity Problem

From the spoke table in CONTEXT.md:

| my_spoke | orient | CW result | CCW result |
|----------|--------|-----------|------------|
| 0 (H)   | H      | NW        | NE         |
| 1 (NE)  | NE     | H         | NW         |
| 2 (NW)  | NW     | NE        | H          |
| 3 (H)   | H      | NW        | NE         |
| 4 (NE)  | NE     | H         | NW         |
| 5 (NW)  | NW     | NE        | H          |

If the rule is "avoid NW," then:
- Spoke 0 (H): CCW gives NE (good), CW gives NW (bad) -- choose CCW
- Spoke 1 (NE): CW gives H (good), CCW gives NW (bad) -- choose CW
- **Spoke 2 (NW): CW gives NE (good), CCW gives H (good) -- BOTH are non-NW!**
- Spoke 3 (H): same as spoke 0 -- choose CCW
- Spoke 4 (NE): same as spoke 1 -- choose CW
- **Spoke 5 (NW): same as spoke 2 -- BOTH are non-NW!**

The ambiguity at NW spokes is real. But Purbhoo says "this restriction rules out one of two possible 120 degree rotations." So either:

**(A)** The rule is not "avoid NW" but something more specific, or
**(B)** NW-oriented edges never actually occur at 3-fold symmetric hexagons during migration.

---

## 4. Re-Reading the Rule: "Some Edge of the Rhombus"

I initially read the rule as being about the shared edge at the pivot vertex. But Purbhoo says "some **edge** of the rhombus," not "the shared edge." A rhombus has four edges. After rotation, the rhombus occupies a new position with edges in two orientations (the rhombus's four edges come in two parallel pairs). The rule might mean: **choose the rotation such that one of the rhombus's edge-pairs becomes horizontal or vertical.**

A 30/150-degree rhombus in the equilateral grid has edges in two of the three orientations {H, NE, NW}. After a 120-degree rotation, both edge orientations shift by one step in the cyclic order H -> NE -> NW -> H (or its reverse). So if the rhombus starts with edges in orientations {X, Y}, it ends with edges in {X+1, Y+1} or {X-1, Y-1} (cyclically).

The three possible starting pairs are: {H, NE}, {NE, NW}, {NW, H}. After rotation:
- {H, NE} -> {NE, NW} or {NW, H}
- {NE, NW} -> {NW, H} or {H, NE}
- {NW, H} -> {H, NE} or {NE, NW}

If the rule is "one resulting pair must include H or NE (not be {NE, NW})":
- From {H, NE}: both results include H or NE -- still ambiguous when starting from certain configurations!

This reframing does not resolve the ambiguity cleanly either. The problem persists for some starting orientations regardless of whether we interpret the rule as being about one edge or all edges.

---

## 5. The Correct Resolution: Proposition 3.1 and the Two-Orientation Constraint

Proposition 3.1 (page 14) provides a critical constraint:

> "All rhombi in [the migrated collection] are in the same nest. Consequently, [the collection] is a transformed diagram."

The proof states: "for single (hexagonal) mosaics this is clear, as the rhombi can only ever be in **two of three orientations** during the migration process, thereby ruling out one of the nests."

This is the key. During migration, a rhombus cycles through only TWO of the three possible orientations, never all three. The direction rule enforces this constraint. For A-to-B migration, the rhombus alternates between orientations that belong to nests A and B, never visiting the orientation associated with nest C.

In the hex_mosaic.py coordinate system, if nest C's orientation is NW, then the direction rule is precisely: **never let the rhombus's edge become NW-oriented.** Since starting from an NW spoke would mean the rhombus is ALREADY in the forbidden orientation, the question is whether this situation arises. If migration is working correctly, the rhombus should never be in a NW-oriented state at a 3-fold symmetric hexagon -- it should only arrive there in H or NE states. The NW spokes (2 and 5) should never be the active spoke during a correctly-operating migration from A to B.

**However**, this depends on correctly identifying which orientation belongs to which nest. The mapping nest-orientation to grid-orientation {H, NE, NW} depends on the hexagon geometry and should be verified empirically from Figure 2.

---

## 6. The Direction Rule's Interaction with Migration Direction

The migration direction determines which nest is the target, which determines which orientation is forbidden, which determines the direction of each 120-degree rotation. The chain is:

1. **Migration target** (e.g., A to B) fixes the compass direction
2. **Compass direction** determines which pair of nest orientations are "allowed" -- specifically, the orientations of the source and target nests
3. **Allowed orientations** determine the excluded orientation
4. At each 3-fold hexagon, **choose the rotation whose resulting edge orientation is in the allowed set**
5. Since the starting orientation is always in the allowed set (by induction on migration steps), one rotation keeps it in the allowed set and one sends it to the excluded orientation. **Exactly one rotation is valid.**

This resolves the NW-spoke ambiguity: if the migration is proceeding correctly, the rhombus never reaches a state where its spoke orientation is in the excluded family. The ambiguity in the table is an artifact of considering all six spokes abstractly, when in practice only four of the six can occur during any given migration.

---

## 7. Summary and Practical Implications

**The direction rule for A-to-B migration is**: at each 3-fold symmetric hexagon, rotate 120 degrees in the direction that keeps the rhombus's edge orientation within the two-element set {orientation_A, orientation_B}, avoiding orientation_C.

**For the implementation**: The fix is NOT a per-spoke lookup table (which would need to handle the impossible NW-spoke case). The fix is:
1. Determine which orientation is forbidden for the given migration (source nest, target nest)
2. At each 3-fold hexagon, compute both rotation results
3. Choose the rotation whose resulting orientation is NOT the forbidden one
4. Assert that exactly one rotation satisfies this (if both do, the rhombus has entered an impossible state -- this is a bug)

**The forbidden orientation** for A-to-B migration needs to be determined from the hexagon geometry. From the hex_mosaic.py spoke mapping and the hexagon structure in Figure 5, this should be one of {H, NE, NW}. My geometric analysis above suggests NW, but this must be verified against Figure 2's worked example -- trace the actual orientation changes of the cyan rhombus through its migration steps and confirm which orientation never appears.

---

## 8. Open Questions for the Vanilla Center Agent

1. Can you verify from Figure 2 which of {H, NE, NW} is the forbidden orientation for the leftward-to-rightward migration shown? The cyan and magenta rhombi's orientations at each step would empirically settle this.

2. The proof of Proposition 3.1 says "two of three orientations" but does not specify which two. For a given migration direction, which two are they? This is determined by the nest geometry but the paper leaves it implicit.

3. In the hex_mosaic.py coordinate system, how do the nest orientations (E_A, N_A), (E_B, N_B), (E_C, N_C) from Figure 5 map to the grid orientations {H, NE, NW}? This mapping is the missing link between Purbhoo's abstract rule and the concrete implementation.
