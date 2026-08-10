#!/usr/bin/env python3
"""Implement Purbhoo's migration bijection: puzzle → mosaic → migrate flock → LR tableau.

The migration operates on the MOSAIC (hexagonal representation), not the puzzle
(triangular representation). The puzzle's 10-labels correspond to "square" tiles
in the mosaic. The FLOCK is formed from corner rhombi at nest A (left boundary),
which are reconstructed from the boundary data.

Reference: Purbhoo, "Puzzles, Tableaux, and Mosaics" (arXiv:0705.1184)
"""

import sys
from pathlib import Path
from copy import deepcopy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "puzzles"))
from lr_coefficients import (
    lr_puzzles, _DELTA_FUNC, _NABLA_BRANCH,
    partition_to_boundary, boundary_to_partition, _bounding_box,
)


def puzzle_to_edge_grid(puzzle):
    """Convert a puzzle to a dict of all edge labels indexed by endpoints.

    Returns dict: {((r1,j1), (r2,j2)): label}
    where (r1,j1) < (r2,j2) lexicographically.
    """
    return puzzle.edge_labels()


def get_horizontal_edges(puzzle):
    """Get all horizontal edge labels, indexed by (level, position)."""
    horiz = {}
    for r, lev in enumerate(puzzle.levels):
        for j, val in enumerate(lev):
            horiz[(r, j)] = val
    return horiz


def trace_rhombus_positions(puzzle):
    """Find all positions where a rhombus tile (KTW 10-piece) sits.

    A rhombus tile in the puzzle spans one △ and one ∇. It's identified
    by the horizontal edge with label 10 between them. Returns a list of
    (level, position) tuples where 10 appears in the levels.
    """
    positions = []
    for r, lev in enumerate(puzzle.levels):
        for j, val in enumerate(lev):
            if val == 10:
                positions.append((r, j))
    return positions


def trace_partition_evolution(puzzle):
    """Track the partition at each even level (pure {0,1} boundary words).

    At odd levels, the 10-label prevents direct partition decoding.
    Returns: list of (level, partition) for even levels only.
    """
    results = []
    for r, lev in enumerate(puzzle.levels):
        if 10 not in lev and len(lev) > 0:
            part = boundary_to_partition(lev)
            results.append((r, part))
        elif len(lev) == 0:
            results.append((r, ()))
    return results


def compute_level_skews(puzzle):
    """For each pair of consecutive even-level partitions, compute the skew shape.

    The boxes removed between consecutive even levels tell us which Young diagram
    boxes are "processed" in that band of the triangle.
    """
    even_parts = trace_partition_evolution(puzzle)
    skews = []
    for i in range(len(even_parts) - 1):
        r1, p1 = even_parts[i]
        r2, p2 = even_parts[i + 1]
        # Boxes in p1 but not in p2
        removed = []
        for row in range(max(len(p1), len(p2))):
            c1 = p1[row] if row < len(p1) else 0
            c2 = p2[row] if row < len(p2) else 0
            for col in range(c2, c1):
                removed.append((row, col))
        skews.append((r1, r2, removed))
    return skews


def direct_bijection_attempt(lam, mu, nu):
    """Attempt to find the puzzle → LR-tableau bijection by tracking
    which boxes of ν are removed at each band of the transfer matrix,
    and comparing with which rows of μ contribute at each band.

    The left boundary bits tell us when we cross a ROW of λ (1-bits).
    The right boundary bits tell us when we cross a ROW of μ (1-bits).
    Each band between consecutive row-crossings processes specific boxes.
    """
    puzzles = lr_puzzles(lam, mu, nu)
    n, m = _bounding_box(lam, mu, nu)
    N = n + m

    lam_bits = partition_to_boundary(lam, n, m)
    mu_bits = partition_to_boundary(mu, n, m)

    # Right boundary is read as right[N-1-r] for row r
    print(f"λ boundary (left):  {lam_bits}")
    print(f"μ boundary (right): {mu_bits}")
    print(f"Right boundary at each row: {[mu_bits[N-1-r] for r in range(N)]}")
    print()

    # The μ rows are crossed when right[N-1-r] = 1
    # This happens at rows where the right boundary has a 1-bit
    mu_row_levels = [r for r in range(N) if mu_bits[N - 1 - r] == 1]
    lam_row_levels = [r for r in range(N) if lam_bits[r] == 1]
    print(f"λ row crossings at levels: {lam_row_levels}")
    print(f"μ row crossings at levels: {mu_row_levels}")
    print()

    # For each puzzle, track which boxes of ν/λ are in which band
    skew_boxes = []
    for i in range(len(nu)):
        lam_i = lam[i] if i < len(lam) else 0
        for j in range(lam_i, nu[i]):
            skew_boxes.append((i, j))

    print(f"Skew boxes ν/λ: {skew_boxes}")
    print()

    # For each puzzle, compute which even-level partition bands
    # contain which skew boxes
    for pi, puzzle in enumerate(puzzles):
        print(f"=== Puzzle {pi + 1} ===")
        skews = compute_level_skews(puzzle)

        # Map each skew box to the band where it's removed
        box_to_band = {}
        for band_idx, (r1, r2, removed) in enumerate(skews):
            for box in removed:
                if box in skew_boxes:
                    box_to_band[box] = (band_idx, r1, r2)

        print(f"  Even-level partition evolution:")
        even_parts = trace_partition_evolution(puzzle)
        for r, p in even_parts:
            print(f"    Level {r}: {p}")

        print(f"  Skew box → band mapping:")
        for box in skew_boxes:
            if box in box_to_band:
                band_idx, r1, r2 = box_to_band[box]
                print(f"    Box {box}: removed in band {band_idx} (levels {r1}→{r2})")

        # Now: the μ row crossings tell us which "content entry" each band corresponds to
        # Band i spans from even_level i to even_level i+1
        # The μ-row crossings within that band determine the entry value
        print(f"  Attempting content assignment:")
        mu_row_counter = 0
        tableau = {}
        for band_idx, (r1, r2, removed) in enumerate(skews):
            # Count μ-row crossings in this band
            crossings = [r for r in mu_row_levels if r1 <= r < r2]
            if crossings:
                mu_row_counter += 1
            for box in removed:
                if box in skew_boxes:
                    tableau[box] = mu_row_counter
        print(f"    Tableau: {tableau}")
        print(f"    Entries at skew boxes: {[tableau.get(b, '?') for b in skew_boxes]}")
        print()


def rhombus_path_analysis(lam, mu, nu):
    """Analyze the 10-positions and their relationship to the bijection.

    Track the R (=10) labels through the puzzle levels and analyze how they
    relate to the LR-tableau entries.
    """
    puzzles = lr_puzzles(lam, mu, nu)

    print(f"Rhombus (10) positions in levels:")
    print()

    for pi, puzzle in enumerate(puzzles):
        print(f"=== Puzzle {pi + 1} ===")
        r_positions = trace_rhombus_positions(puzzle)
        print(f"  R positions: {r_positions}")

        # For each R at an odd level, compute the "carried box"
        # by expanding 10 → (1,0) vs (0,1) and comparing partitions
        for r, j in r_positions:
            lev = puzzle.levels[r]
            upper = tuple(x for v in lev for x in ([1, 0] if v == 10 else [v]))
            lower = tuple(x for v in lev for x in ([0, 1] if v == 10 else [v]))
            u_part = boundary_to_partition(upper)
            l_part = boundary_to_partition(lower)

            # Find the carried box
            for row in range(max(len(u_part), len(l_part))):
                uv = u_part[row] if row < len(u_part) else 0
                lv = l_part[row] if row < len(l_part) else 0
                if uv > lv:
                    print(f"  R at level {r}, pos {j}: carries box ({row},{lv})"
                          f"  [upper={u_part}, lower={l_part}]")
        print()


if __name__ == "__main__":
    lam, mu, nu = (2, 1), (2, 1), (3, 2, 1)

    print("=" * 60)
    print(f"Puzzle → LR-Tableau Bijection Analysis")
    print(f"c^{nu}_{{{lam},{mu}}} = {len(lr_puzzles(lam, mu, nu))}")
    print("=" * 60)
    print()

    print("LR Tableaux (computed by hand):")
    print("  T1: (0,2)=1, (1,1)=1, (2,0)=2  [reading word: 1,1,2]")
    print("  T2: (0,2)=1, (1,1)=2, (2,0)=1  [reading word: 1,2,1]")
    print()

    print("-" * 60)
    print("ANALYSIS 1: Rhombus path tracking")
    print("-" * 60)
    rhombus_path_analysis(lam, mu, nu)

    print("-" * 60)
    print("ANALYSIS 2: Direct bijection via partition band decomposition")
    print("-" * 60)
    direct_bijection_attempt(lam, mu, nu)

    print("-" * 60)
    print("ANALYSIS 3: Level-by-level comparison")
    print("-" * 60)
    puzzles = lr_puzzles(lam, mu, nu)
    for pi, puzzle in enumerate(puzzles):
        print(f"\nPuzzle {pi + 1}:")
        print(puzzle)
