# Orientation Cycle Analysis: CW vs CCW at Each Spoke

## 1. The Core Structure

The spoke-to-orientation map creates a repeating pattern with period 2:

```
Spoke:   0    1    2    3    4    5
Orient:  H    NE   NW   H    NE   NW
```

CW rotation adds 2 to the spoke index (mod 6); CCW subtracts 2. This means every rotation moves the rhombus two spokes around the hexagon. The resulting orientation table, reproduced from CONTEXT.md:

| spoke | orient | CW -> spoke (orient) | CCW -> spoke (orient) |
|-------|--------|----------------------|----------------------|
| 0     | H      | 2 (NW)               | 4 (NE)              |
| 1     | NE     | 3 (H)                | 5 (NW)              |
| 2     | NW     | 4 (NE)               | 0 (H)               |
| 3     | H      | 5 (NW)               | 1 (NE)              |
| 4     | NE     | 0 (H)                | 2 (NW)              |
| 5     | NW     | 1 (NE)               | 3 (H)               |

The structure has Z/3 symmetry: spokes {0,1,2} and {3,4,5} behave identically under the orient map, just shifted by 3.


## 2. Purbhoo's Direction Rule: What It Actually Says

Section 3.1, page 13: "If we see (b), the hexagon has a 3-fold rotational symmetry, so we rotate the tiling 120 degrees, *in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.* (This restriction rules out one of the two possible 120 degree rotations.)"

The critical terms are "horizontal" and "vertical." Purbhoo warns on page 4 that all directions are relative: "When we speak of 'up', 'down', 'left' and 'right', these are defined relative to a chosen 'right'." And on page 13: "In describing migration we will assume that our diagrams are turned so that the direction of migration is pointing directly to our right."

For A-to-B migration, the migration direction is "to the right." In Purbhoo's mosaic geometry, the three tile shapes are:
- (a) equilateral triangle, side 1
- (b) unit square (this is the H-rhombus in the grid: angles 90/90)
- (c) 30/150-degree rhombus (this becomes NW or NE in the grid)

Now, Purbhoo's mosaic sits inside the hexagon A'AB'BC'C with 150-degree angles at A, B, C. The three nests have orientations defined by the hexagon's edge directions. After the "straightening" operation (page 11) that maps mosaic to equilateral-triangle puzzle, the squares become rhombi and the directions rotate.

**The key geometric fact**: In the equilateral triangle grid used by `hex_mosaic.py`, the three edge orientations are:
- **H** (horizontal): edges like ((r,j), (r,j+1)), running left-right
- **NE** (northeast): edges like ((r,j), (r+1,j)), running lower-left to upper-right
- **NW** (northwest): edges like ((r,j), (r+1,j-1)), running lower-right to upper-left

When Purbhoo says "horizontal or vertical" with migration direction pointing right, he means edges that are aligned with or perpendicular to the migration direction. In the equilateral grid with rightward migration (increasing j):
- **H edges run parallel to the migration direction** -- these are "horizontal"
- **NE edges run at 60 degrees** to horizontal -- in Purbhoo's coordinate system with the page rotated, these are "vertical" (perpendicular to migration)
- **NW edges run at 120 degrees** to horizontal -- these are neither horizontal nor vertical; they are the diagonal direction

Therefore: **after a 120-degree rotation, the resulting spoke orientation must be H or NE, not NW.**


## 3. Applying the Rule to Each Spoke

### Unambiguous cases (4 of 6):

**Spoke 0 (H):** CW gives NW (forbidden), CCW gives NE (allowed). **Answer: CCW (shift=2).**

**Spoke 1 (NE):** CW gives H (allowed), CCW gives NW (forbidden). **Answer: CW (shift=4).**

**Spoke 3 (H):** CW gives NW (forbidden), CCW gives NE (allowed). **Answer: CCW (shift=2).**

**Spoke 4 (NE):** CW gives H (allowed), CCW gives NW (forbidden). **Answer: CW (shift=4).**

These four are cleanly determined. The pattern: H-spokes go CCW, NE-spokes go CW.

### Ambiguous cases (2 of 6):

**Spoke 2 (NW):** CW gives NE (allowed), CCW gives H (allowed). **Both results are non-NW.**

**Spoke 5 (NW):** CW gives NE (allowed), CCW gives H (allowed). **Both results are non-NW.**

This is the NW ambiguity. The direction rule as stated ("end up horizontal or vertical") does not distinguish between the two options when starting from NW, because both CW and CCW produce acceptable outcomes.


## 4. Resolving the NW Ambiguity

There are four approaches to resolving this. I will examine each.

### Approach A: The "rightward progress" principle

Purbhoo says the rhombus "moves roughly to the right (although sometimes it moves strictly vertically)." For rightward progress, we want the rotation that moves the rhombus position to higher j (or at least does not decrease j). This is what the old `direction_score` heuristic tried to do -- but it failed because the positional effect depends on the local tiling configuration, not just the spoke index.

However, we can examine the *typical* positional effect. Starting from spoke 2 (NW), with pivot vertex at (r, j+1):
- CW to spoke 4 (NE): new position via `_spoke_to_position(r, j+1, 4, 'ne')` = (r-1, j). The j-coordinate **decreased** by 1.
- CCW to spoke 0 (H): new position via `_spoke_to_position(r, j+1, 0, 'h')` = (r-1, j+1). The j-coordinate is **unchanged**.

For spoke 5 (NW), pivot vertex at ... but wait -- spoke 5 is never reached from `_pivot_vertex`. The pivot function only returns spokes 2, 3, or 4. Let me verify this.

From `_pivot_vertex`:
- orient='nw' returns spoke 2
- orient='h' returns spoke 3
- orient='ne' returns spoke 4

**Spokes 0, 1, 3, 5 are never the starting spoke for a pivot.** Wait -- spoke 3 IS returned for orient='h'. Let me re-examine.

Actually, checking the `_pivot_vertex` function more carefully:
- orient='nw': returns center (r, j+1), spoke **2**
- orient='h': returns center (r+1, j+1), spoke **3**
- orient='ne': returns center (r+1, j+1), spoke **4**

So the starting spokes are always 2, 3, or 4. This means spokes 0, 1, and 5 are never initial states -- they can only be reached as *results* of rotations. The ambiguity therefore only applies to **spoke 2**, not to spoke 5.

This is a significant structural observation. The `_MIGRATION_SHIFT` table has entries for all 6 spokes, but `_pivot_vertex` only ever produces spokes 2, 3, 4. The entries for spokes 0, 1, 5 in `_MIGRATION_SHIFT` are dead code.

### Approach B: Tracing Figure 2

Figure 2 of Purbhoo shows a complete migration from left to right with three rhombi (labels 1, 1, 2). The figure has 14 subfigures arranged in 4 rows. Each subfigure shows one rotation step. By tracking the orientation of the migrating rhombus (the colored one) through each step, we could empirically determine which direction is chosen when starting from NW.

From the figure, the pink and cyan rhombi undergo several 120-degree rotations as they move rightward. The shapes of the rhombi change between frames, which corresponds to orientation changes. Unfortunately, the printed resolution makes it difficult to distinguish NW from NE orientations precisely in every frame. But the overall trajectory is consistently rightward, with each step either maintaining or increasing the j-coordinate. This favors **CCW for spoke 2** (which preserves j) over CW (which decreases j).

### Approach C: The orientation cycle interpretation

The HANDOFF.md suggests a cycle NW -> H -> NE -> NW that advances "rightward." If this cycle governs the direction rule, then:
- NW should rotate to H (CCW for spoke 2)
- H should rotate to NE (CCW for spoke 3)
- NE should rotate to ... NW? But the rule forbids NW!

This cycle interpretation breaks down because it requires producing NW from NE, which contradicts the "no NW result" rule. The cycle interpretation from HANDOFF.md is **inconsistent** with Purbhoo's stated rule.

### Approach D: Algebraic consistency with _MIGRATION_SHIFT

The existing `_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}` encodes:
- shift=2 means CCW
- shift=4 means CW

So the table says: spoke 2 -> CCW (to spoke 0, orient H). This is consistent with Approach A (rightward progress favoring j preservation) and with the 4-spoke analysis above (H and NW use CCW, NE uses CW).


## 5. The Definitive Rule

For A-to-B migration (rightward), the direction at each spoke is:

| spoke | orient | direction | result spoke | result orient |
|-------|--------|-----------|-------------|---------------|
| 2     | NW     | **CCW**   | 0           | H             |
| 3     | H      | **CCW**   | 1           | NE            |
| 4     | NE     | **CW**    | 0           | H             |

Since `_pivot_vertex` only ever produces spokes 2, 3, 4, the rule simplifies to:

```
direction = CW  if my_spoke == 4  (NE)
            CCW if my_spoke in {2, 3}  (NW or H)
```

Or equivalently: **CW only for NE spokes; CCW for everything else.**

This matches the current code at line 979:
```python
direction = 1 if my_spoke in (1, 4) else -1
```

Since spoke 1 is never produced by `_pivot_vertex`, including it is harmless. The code is correct for the three spokes that actually matter (2, 3, 4).


## 6. Does the Direction Depend on Migration Direction (A->B vs B->A)?

Yes, it must. Purbhoo's rule is relative to the migration direction. For B-to-A migration (leftward), the roles of "horizontal" and "vertical" remain the same (they are defined by the page orientation), but "rightward progress" reverses. The pivot vertex would be chosen differently (the vertex with *smaller* j), and the resulting spoke set would change.

Theorem 2 states: "The inverse operation of migration of a flock F from nest A to nest B is the migration of F' from B to A." The direction rule for B-to-A would involve choosing the vertex that leads leftward (decreasing j), and the spoke indices produced by a hypothetical `_pivot_vertex_left` function would differ.

For the current implementation, which only handles A-to-B, the rule as stated in Section 5 is sufficient. If B-to-A is needed later, the pivot vertex selection and direction table would need to be derived separately.


## 7. Is `_MIGRATION_SHIFT` Correct?

The table `_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}` is **correct for the three spokes that matter** (2, 3, 4). For the three spokes that `_pivot_vertex` never produces (0, 1, 5), the entries are best-effort guesses but irrelevant in practice.

However, the current `migrate_single` code at line 979 does NOT use `_MIGRATION_SHIFT`. It uses:
```python
direction = 1 if my_spoke in (1, 4) else -1
```

This is equivalent to `_MIGRATION_SHIFT` for the relevant spokes:
- spoke 2: direction=-1 (CCW), equivalent to shift=2
- spoke 3: direction=-1 (CCW), equivalent to shift=2
- spoke 4: direction=+1 (CW), equivalent to shift=4

So the direction logic in `migrate_single` is already consistent with `_MIGRATION_SHIFT` and with Purbhoo's rule. **The direction rule is not the bug.**


## 8. Where the Bug Must Be

If the direction rule is correct, the collisions at rhombi 3-6 must come from elsewhere. Candidates:
1. **Domain setup**: Missing or extra edges near nest B causing incorrect 2-fold/3-fold classification
2. **Termination condition**: `_is_stuck` or `_is_in_nest_b` failing to detect the correct stopping point
3. **Pivot vertex selection**: `_pivot_vertex` choosing the wrong leading vertex in edge cases
4. **ensure_hex_cell**: Filling edges with label 1 that should be label 10 (or vice versa), corrupting the local tiling state

The direction rule analysis is clean. The NW ambiguity resolves in favor of CCW (the code's current choice). The remaining bug is structural, not directional.


## 9. Note on Vanilla Findings

The vanilla center findings file (`vanilla_findings.md`) does not exist yet. This analysis was conducted independently. When the vanilla analysis becomes available, a key point to cross-check is whether it concurs that the direction rule is correct and the bug lies in domain setup or termination, not in CW/CCW selection.
