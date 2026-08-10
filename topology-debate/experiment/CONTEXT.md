# Shared Context: Purbhoo Migration Direction Rule Problem

## The Goal
Implement Purbhoo's migration bijection (arXiv:0705.1184, Section 3.1) correctly.
The bijection goes: KTW puzzle -> mosaic -> flock on nest A -> migrate A to B -> read LR-tableau at B.
This is the ONLY known direct bijection from puzzles to LR tableaux.

## The Bug
File: `~/git/puzzles/hex_mosaic.py`, function `migrate_single` (line ~915)

**Symptom**: When migrating 6 rhombi from nest A to nest B, rhombi 3-6 sometimes collide (land at the same position as an earlier rhombus). Expected: all 6 land at 6 distinct nest-B positions.

**Root cause**: The 120-degree rotation direction choice uses a POSITION-BASED heuristic:
```python
def direction_score(p):
    pr, pj = p
    return (1 if pr >= 0 else 0, pj, pr)
```
This picks whichever rotation moves the rhombus more to the right (higher j). When both options have the same j, it makes an arbitrary tiebreak that's sometimes wrong.

## What Purbhoo's Paper Actually Says (page 13)

> "If we see (b), the hexagon has a 3-fold rotational symmetry, so we rotate the tiling 120 degrees, **in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.** (This restriction rules out one of two the possible 120 degree rotations.)"

This is THE rule. But what does "horizontal" and "vertical" mean in Purbhoo's coordinate system?

## Purbhoo's Coordinate System

The hexagon A'AB'BC'C has:
- Angles 150 degrees at A, B, C and 90 degrees at A', B', C'
- Side lengths: A'A = B'B = C'C = d and AB' = BC' = CA' = n-d
- Three nests at corners A, B, C
- Orientations: (E_A, N_A), (E_B, N_B), (E_C, N_C) defined by the hexagon edge directions

Mosaic tiles:
- (a) equilateral triangle, side 1
- (b) square, side 1  
- (c) rhombus, side 1, angles 30 and 150 degrees

**CRITICAL**: Purbhoo warns (page 4): "all directions in this paper are used in a relative, rather than in an absolute sense. When we speak of 'up', 'down', 'left' and 'right', these are defined relative to a chosen 'right'; hence the reader may sometimes need to rotate the page."

For migration from A to B, the direction is "to the right" and diagrams should be rotated so the migration direction points right.

## The hex_mosaic.py Coordinate System

Uses an equilateral triangle grid with:
- Face coordinates: (type, r, j) where type is 'up' or 'down'
- Edge coordinates: pairs of vertices ((r1,j1), (r2,j2))
- Three edge orientations: H (horizontal), NE (northeast), NW (northwest)
- Vertices at hex centers have 6 "spokes" (interior edges), numbered 0-5

Spoke-to-orientation mapping:
```
_spoke_to_orient: 0=H, 1=NE, 2=NW, 3=H, 4=NE, 5=NW
```

Rotation mechanics:
- CW 120 degrees: my_spoke -> (my_spoke + 2) % 6
- CCW 120 degrees: my_spoke -> (my_spoke - 2) % 6
- 180 degrees: my_spoke -> (my_spoke + 3) % 6

## The Direction Rule Analysis (Partially Done)

For A->B migration, starting from spoke `my_spoke`, the two rotation options give:

| my_spoke | orient | CW -> spoke | CW orient | CCW -> spoke | CCW orient |
|----------|--------|-------------|-----------|--------------|------------|
| 0        | H      | 2           | NW        | 4            | NE         |
| 1        | NE     | 3           | H         | 5            | NW         |
| 2        | NW     | 4           | NE        | 0            | H          |
| 3        | H      | 5           | NW        | 1            | NE         |
| 4        | NE     | 0           | H         | 2            | NW         |
| 5        | NW     | 1           | NE        | 3            | H          |

If Purbhoo's rule means "the resulting orientation should be H or NE (not NW)":
- spoke 0 (H): CCW gives NE (good) vs CW gives NW (bad) -> CCW
- spoke 1 (NE): CW gives H (good) vs CCW gives NW (bad) -> CW
- spoke 2 (NW): CW gives NE vs CCW gives H -> BOTH are "good"! Ambiguous!
- spoke 3 (H): CCW gives NE (good) vs CW gives NW (bad) -> CCW
- spoke 4 (NE): CW gives H (good) vs CCW gives NW (bad) -> CW
- spoke 5 (NW): CW gives NE vs CCW gives H -> BOTH are "good"! Ambiguous!

**The NW spokes (2 and 5) are ambiguous** under this interpretation. This is where the bug lives.

## Existing Code Workaround

The code has a static table that was meant to resolve this:
```python
_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}
```
But `migrate_single` doesn't use this table! It uses the `direction_score` heuristic instead.

## Key Questions to Resolve

1. What does "horizontal or vertical" mean precisely in Purbhoo's mosaic coordinates for A->B migration?
2. How do the three rhombus types (30/150 degree) in the mosaic map to NW/H/NE in the equilateral grid?
3. For NW-oriented spokes, which 120-degree direction should be chosen?
4. Is the direction rule even about the RESULTING orientation, or about which EDGE of the rhombus becomes horizontal/vertical?
5. How should the nest-B domain be set up before migration begins?
6. What is the correct termination condition?

## Figure 2 (Worked Example)

Figure 2 shows a complete migration of a flock from left to right. The flock has labels 1, 1, 2 (three rhombi). Each row of subfigures shows the rhombus moving through several rotation steps. The migration direction is to the right.

Tracing the orientation changes visible in the figure would empirically determine the direction rule.

## Test Case

Boundaries: lambda = mu = (2,1), nu = (3,2,1), N = 6
Expected: 6 rhombi migrate to 6 distinct positions [(1,4),(1,5),(1,6),(3,2),(3,3),(5,0)]
Current: Rhombus 1 works correctly, rhombi 3-6 collide.

## Key References in this Directory

- `HANDOFF.md` — The bug description and reproduction steps
- `FINDINGS.md` — Hypotheses S1 and S5 about migration
- `RESULTS.md` — Ring vs star topology experiment results

## Files in ~/git/puzzles/

- `hex_mosaic.py` — The migration implementation (the code with the bug)
- `mosaic.py` — Base mosaic class with edge-label storage
- `lr_coefficients.py` — Generates KTW puzzles
- `seed-papers/purbhoo-2007-puzzles-tableaux-mosaics.pdf` — The paper
