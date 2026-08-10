# Agent 1 (Paper Exegete): Purbhoo's Direction Rule -- A Close Reading

## Summary of Finding

Purbhoo's "horizontal or vertical" direction rule for 120-degree rotations is **not about the grid orientation labels (H, NE, NW) in the equilateral triangle grid**. It refers to directions in the **mosaic's physical geometry**, which is a hexagon embedded in the plane with specific edge directions determined by the nest orientations. The rule must be interpreted relative to the **migration direction**, which Purbhoo explicitly tells us to treat as "rightward." The current implementation's failure stems from conflating grid-level orientation names with the mosaic's physical directions.

---

## 1. The Warning on Page 4: Directions Are Relative

Purbhoo writes (page 4, emphasis mine):

> "We would also like to warn the reader that all directions in this paper are used in a **relative**, rather than in an absolute sense. When we speak of 'up', 'down', 'left' and 'right', these are defined relative to a chosen 'right'; hence the reader may sometimes need to rotate the page."

And further:

> "Compass directions 'north', 'east', 'south' and 'west', are used even more loosely. These are defined relative to a chosen basis for R^2, which may not be orthogonal or even have the expected orientation."

This is not a throwaway remark. It is the key to the entire direction rule. The mosaic hexagon A'AB'BC'C has three pairs of parallel sides, and the "directions" E, N, -E, -N are defined **by the orientation of the nest**, not by any absolute compass. The unit vectors E_A, N_A, E_B, N_B, E_C, N_C point along the hexagon's edges:

- E_A and N_A point along edges AA' and AB' respectively
- E_B and N_B point along edges BB' and BC' respectively  
- E_C and N_C point along edges CC' and CA' respectively

These are at 150 degrees to each other at the nest corners (A, B, C) and at 90 degrees at the primed corners (A', B', C').

## 2. What "Horizontal" and "Vertical" Mean

On page 13, Purbhoo states:

> "If we see (b), the hexagon has a 3-fold rotational symmetry, so we rotate the tiling 120 degrees, **in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.** (This restriction rules out one of two the possible 120 degree rotations.)"

The critical phrase is "exactly horizontal, or exactly vertical." This refers to directions in the **diagram as drawn**, where the migration direction is oriented to the right. Purbhoo instructs us on page 13:

> "In describing migration we will assume that our diagrams are turned so that the direction of migration is pointing directly to our right."

So "horizontal" means **parallel to the migration direction**, and "vertical" means **perpendicular to the migration direction**.

### The Three Edge Directions of a Mosaic

A mosaic has three edge directions, corresponding to the three pairs of parallel sides of the hexagonal domain. In the equilateral triangle grid, these correspond to:

- **H (horizontal)**: edges connecting vertices at the same row level
- **NE (northeast)**: edges going up-right  
- **NW (northwest)**: edges going up-left

But in the mosaic's physical embedding, these three directions are **not** at 60 degrees to each other as in a regular triangular grid. The mosaic hexagon has angles of 150 degrees at A, B, C and 90 degrees at A', B', C'. The actual directions of the mosaic edges depend on which nest we are near and which migration direction we are considering.

### For A-to-B Migration

The migration direction from A to B is "south" in Purbhoo's notation (see Corollary 3.2, page 14: "Let F migrate from A to B (south)"). When we rotate the diagram so this direction points right, the key question is: **which of the three edge directions become horizontal and which become vertical?**

The hexagon's geometry determines this. The migration direction from A to B is along the edge from nest A toward the AB' side. In the diagram rotated so migration points right:

- One edge direction aligns **horizontally** (parallel to migration direction)
- One edge direction aligns **vertically** (perpendicular to migration direction)
- The third edge direction is at an **oblique angle** (neither horizontal nor vertical)

The direction rule then says: after the 120-degree rotation, the rhombus's shared edge should end up in one of the first two directions (horizontal or vertical), **not** the oblique third direction.

## 3. Mapping to Grid Orientations: The Nest-Dependent Correspondence

Here is the crucial insight that the current code misses. The mapping from {H, NE, NW} (grid labels) to {horizontal, vertical, oblique} (mosaic geometry) **depends on which nest the migration goes from and to**.

For A-to-B migration (south, oriented rightward in the diagram):
- The sides of the hexagon near the A-B corridor determine which grid directions align with horizontal and vertical.
- From Figure 5 (page 10), the hexagon A'AB'BC'C has: side AB' running from A to B' with length n-d, and sides AA' (length d) and B'B (length d) flanking the nests.

Looking at Figure 2 (page 5), which shows a complete migration from left (nest A) to right (nest B), we can trace the rhombus orientations at each step:

**Row 1 of Figure 2** (the "2" rhombus, colored cyan):
- Starts at nest A in a tilted orientation
- After first move: changes orientation
- After second move: changes again
- After third move: changes again (reaches nest B)

**Row 2** (the "1" rhombus below, colored magenta):
- Follows a similar sequence of orientation changes

**Row 3** (the bottom "1" rhombus):
- Two steps to reach nest B

**Row 4** (final state with both rhombi settled at nest B)

The orientations visible in Figure 2 show the rhombi cycling through directions as they move rightward. At each 3-fold rotation step, the chosen direction is the one that makes the rhombus's edge horizontal or vertical **in the diagram as drawn** (with migration pointing right).

## 4. The Ambiguity at NW Spokes: Why Both Options Look "Good"

From the CONTEXT.md analysis, when the current spoke is NW (spokes 2 and 5):

| my_spoke | orient | CW -> spoke | CW orient | CCW -> spoke | CCW orient |
|----------|--------|-------------|-----------|--------------|------------|
| 2        | NW     | 4           | NE        | 0            | H          |
| 5        | NW     | 1           | NE        | 3            | H          |

Both CW and CCW produce a "non-NW" result (NE or H). This is precisely because NW is the **oblique** direction -- the one that is neither horizontal nor vertical. A rhombus starting with an NW-oriented shared edge has that edge in the "forbidden" direction, so EITHER rotation takes it to an allowed direction.

But the rule is not "avoid NW." The rule is "end up horizontal OR vertical." These are two different targets, and Purbhoo says the rule "rules out one of two possible 120-degree rotations." This means for every starting orientation, exactly one of CW/CCW is correct. Therefore, the NW case must also have a unique answer.

### Resolution: "Horizontal" vs "Vertical" Are Not Interchangeable

The resolution is that "horizontal" and "vertical" are **specific directions**, and for a given starting spoke, only ONE of the two allowed targets is reachable by each rotation. The rule picks the rotation that hits EITHER target -- but since the two rotations from NW go to H and NE respectively, the question becomes: **is H the horizontal direction or the vertical direction?** And **is NE the horizontal or vertical?**

For A-to-B migration with migration pointing right:
- If H (grid-horizontal) aligns with the physical horizontal (migration direction), then H is "horizontal" and NE is "vertical"
- The NW direction is the oblique one at roughly 120 degrees to horizontal

Under this interpretation:
- **From NW spoke (2 or 5)**: CCW gives H ("horizontal"), CW gives NE ("vertical"). Both satisfy the rule! 

This seems to reproduce the ambiguity. But wait -- the rule says "some edge of the rhombus." A rhombus has four edges. The shared edge (the spoke) changes orientation, but the **other edges** of the rhombus also have orientations. The rule might be about ANY edge of the rhombus, not just the shared/spoke edge.

## 5. Re-reading the Rule: "Some Edge of the Rhombus"

The exact text is: "in whichever way causes **some edge of the rhombus** to end up exactly horizontal, or exactly vertical."

A rhombus (30/150 degrees) has four edges in two parallel pairs. In the mosaic, the rhombus's edges are in TWO of the three directions (not all three). The "shared edge" (the spoke) is one direction; the other two edges are in a second direction.

After a 120-degree rotation, the spoke moves to a new direction, and the other pair of edges also moves. The rule says that after rotation, at least one edge (from either pair) should be horizontal or vertical.

Since the rhombus has edges in exactly two of three directions, and after rotation these map to two OTHER directions, we need: at least one of these two resulting directions to be horizontal or vertical.

Given three directions {horizontal, vertical, oblique}, and the rhombus occupies two of them after rotation:
- If neither resulting direction is horizontal or vertical, both are oblique -- impossible since there's only one oblique direction
- So one of the two resulting directions is always horizontal or vertical

This would mean the rule is **always satisfied** by both rotations, which contradicts Purbhoo's claim that it "rules out one of two possible 120-degree rotations."

## 6. The Correct Interpretation: The Rule is About the SPOKE Edge Specifically

The contradiction in Section 5 suggests that "some edge" does not literally mean any of the four edges, but specifically the **shared edge** between the rhombus and the hexagon center -- i.e., the spoke edge. This is the edge that determines which face the rhombus "belongs to" after the rotation.

Under this reading: the spoke edge must become horizontal or vertical. Since:
- Starting from H (spokes 0, 3): CW gives NW, CCW gives NE. Only CCW (-> NE) gives vertical. **Choose CCW.**
- Starting from NE (spokes 1, 4): CW gives H, CCW gives NW. Only CW (-> H) gives horizontal. **Choose CW.**
- Starting from NW (spokes 2, 5): CW gives NE, CCW gives H. CW gives vertical, CCW gives horizontal. **Both work.**

The NW ambiguity persists. This means we need additional information to resolve the NW case.

## 7. The Migration Direction Breaks the Tie

Purbhoo says (page 13-14): "In this way, the rhombus moves roughly to the right (although sometimes it moves strictly vertically)."

This gives us additional constraint: the rhombus should move **roughly to the right** after each step. The parenthetical "(although sometimes it moves strictly vertically)" tells us that vertical movement is acceptable but horizontal movement is the default.

For NW spokes, the two options are:
- CW -> NE (vertical in the diagram): the rhombus moves vertically (up or down)
- CCW -> H (horizontal in the diagram): the rhombus moves horizontally (to the right)

Since horizontal (rightward) is the preferred direction of travel, the default for NW should be **CCW -> H** (move rightward). The vertical option (CW -> NE) would be the exception, occurring only "sometimes."

But this interpretation would make the `_MIGRATION_SHIFT` table:
```
{0: 2(CCW), 1: 4(CW), 2: 2(CCW), 3: 2(CCW), 4: 4(CW), 5: 2(CCW)}
```
Which is exactly the table in the current code! The table `_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}` (where shift=2 is CCW and shift=4 is CW) is consistent with this interpretation:
- H spokes (0, 3): CCW -> NE (vertical). Correct.
- NE spokes (1, 4): CW -> H (horizontal). Correct.
- NW spokes (2, 5): CCW -> H (horizontal). This is the "move rightward" default.

**The `_MIGRATION_SHIFT` table appears to be correct.** The bug is that `migrate_single` does not USE this table -- it uses the `direction_score` heuristic instead.

## 8. Actionable Conclusion

The fix is straightforward: **replace the `direction_score` heuristic with the `_MIGRATION_SHIFT` lookup table** that already exists in the code but is unused.

Specifically, in `migrate_single`, lines 967-983 should become:

```python
# 3-fold: use Purbhoo's direction rule (H/V constraint)
shift = _MIGRATION_SHIFT[my_spoke]
if shift == 4:  # CW
    hex_rotate(hm, cr, cj, shift=4)
    new_spoke = (my_spoke + 2) % 6
else:  # CCW (shift == 2)
    hex_rotate(hm, cr, cj, shift=2)
    new_spoke = (my_spoke - 2) % 6
```

However, there is an important caveat: **the NW case might not always choose CCW.** The rule might need to be context-dependent (choosing CW -> NE "sometimes" when the geometry forces vertical movement). The parenthetical "although sometimes it moves strictly vertically" suggests that the vertical option is needed in some configurations -- perhaps when the rightward path is blocked.

### Open question for other agents

Does the NW direction rule depend on **where the rhombus is** in the mosaic (near the top/bottom boundary vs. interior), or is it always CCW? Tracing Figure 2's worked example carefully could empirically resolve this. The `_MIGRATION_SHIFT` static table assumes it is ALWAYS CCW for NW spokes, which may be wrong for boundary cases.

## 9. Summary of Grid-to-Mosaic Direction Mapping (A->B Migration)

| Grid Label | Mosaic Direction (diagram rotated so migration = right) | Purbhoo Category |
|-----------|--------------------------------------------------------|-----------------|
| H         | Horizontal (parallel to migration)                      | "exactly horizontal" |
| NE        | Vertical (perpendicular to migration)                   | "exactly vertical" |
| NW        | Oblique (~120 degrees from horizontal)                  | Neither -- forbidden as target |

Under this mapping:
- From H spoke: must go to NE (vertical). Use CCW.
- From NE spoke: must go to H (horizontal). Use CW.
- From NW spoke: can go to either H (CCW) or NE (CW). Default is CCW (rightward motion), but CW (vertical motion) may be needed at boundaries.

The `direction_score` heuristic fails because it tries to select the rotation that moves the rhombus rightward by POSITION, which is an approximation that breaks down on ties. The correct approach is to select by RESULTING ORIENTATION, using the `_MIGRATION_SHIFT` table.
