"""
Game of Life as co-Kleisli composition on the Grid comonad.

The Grid comonad wraps a 2D toroidal grid with a focused position.
Co-Kleisli arrows are functions  Grid[A] -> B  that read from the
ambient neighborhood to produce a single value.  Composing them
(co-Kleisli composition) runs the first arrow *everywhere*, building
a new ambient grid for the second arrow to read from.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, TypeVar

A = TypeVar("A")
B = TypeVar("B")
C = TypeVar("C")


@dataclass
class Grid:
    """A toroidal 2D grid with a focused position."""

    cells: list[list[int]]
    focus: tuple[int, int]

    @property
    def height(self) -> int:
        return len(self.cells)

    @property
    def width(self) -> int:
        return len(self.cells[0])

    # ── Comonad operations ──────────────────────────────

    def extract(self) -> int:
        """Return the value at the focused position."""
        r, c = self.focus
        return self.cells[r][c]

    def duplicate(self) -> Grid:
        """
        Grid[A] -> Grid[Grid[A]]

        At every position, place the entire grid focused there.
        We represent Grid[Grid[A]] as a Grid whose 'cells' are
        Grid objects — but for efficiency, we just refocus the
        same underlying data (toroidal, so it's safe).
        """
        rows = []
        for r in range(self.height):
            row = []
            for c in range(self.width):
                row.append(Grid(self.cells, (r, c)))
            rows.append(row)
        return Grid(rows, self.focus)

    def extend(self, f: Callable[[Grid], any]) -> Grid:
        """
        (Grid[A] -> B) -> Grid[A] -> Grid[B]

        Apply f at every position, each time with the grid
        focused there. This is  fmap f . duplicate.
        """
        rows = []
        for r in range(self.height):
            row = []
            for c in range(self.width):
                refocused = Grid(self.cells, (r, c))
                row.append(f(refocused))
            rows.append(row)
        return Grid(rows, self.focus)

    # ── Helpers ─────────────────────────────────────────

    def peek(self, dr: int, dc: int) -> int:
        """Read a value at an offset from the focus (toroidal)."""
        r, c = self.focus
        return self.cells[(r + dr) % self.height][(c + dc) % self.width]

    def __repr__(self) -> str:
        if not self.cells or not isinstance(self.cells[0][0], int):
            return f"Grid({self.height}x{self.width}, focus={self.focus})"
        lines = []
        for r, row in enumerate(self.cells):
            parts = []
            for c, val in enumerate(row):
                marker = f"[{val}]" if (r, c) == self.focus else f" {val} "
                parts.append(marker)
            lines.append("".join(parts))
        return "\n".join(lines)


# ── Co-Kleisli arrows ──────────────────────────────────

NEIGHBORS = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, -1),           (0, 1),
    (1, -1),  (1, 0),  (1, 1),
]


def count_neighbors(grid: Grid) -> int:
    """Grid[Cell] -> Int: count living neighbors of the focused cell."""
    return sum(grid.peek(dr, dc) for dr, dc in NEIGHBORS)


def life_rule(grid: Grid) -> int:
    """
    Grid[Cell] -> Cell: Conway's Game of Life rule.

    A cell is born or survives based on its neighbor count:
      - A living cell with 2 or 3 neighbors survives
      - A dead cell with exactly 3 neighbors is born
      - All other cells die (loneliness or overcrowding)
    """
    me = grid.extract()
    n = count_neighbors(grid)
    if me == 1 and n in (2, 3):
        return 1  # survival
    if me == 0 and n == 3:
        return 1  # birth
    return 0  # death




# ── Co-Kleisli composition ─────────────────────────────

def co_kleisli_compose(
    f: Callable[[Grid], any],
    g: Callable[[Grid], any],
) -> Callable[[Grid], any]:
    """
    Compose two co-Kleisli arrows:

        g <=< f  =  g . extend f

    First run f at every position (building a new ambient grid),
    then apply g which reads from that new environment.
    """
    def composed(grid: Grid) -> any:
        extended = grid.extend(f)  # f everywhere -> Grid[B]
        return g(extended)         # g reads from Grid[B]
    return composed


# ── Simulation ──────────────────────────────────────────

def step(grid: Grid) -> Grid:
    """Advance one generation by extending the life rule everywhere."""
    return grid.extend(life_rule)


def run(grid: Grid, generations: int) -> None:
    """Run the Game of Life, printing each generation."""
    for gen in range(generations + 1):
        print(f"Generation {gen}:")
        print(grid)
        print()
        if gen < generations:
            grid = step(grid)


# ── Demo ────────────────────────────────────────────────

def make_grid(pattern: list[str], focus: tuple[int, int] = (0, 0)) -> Grid:
    """Create a grid from a list of strings ('.' = dead, '#' = alive)."""
    cells = []
    for row in pattern:
        cells.append([1 if ch == "#" else 0 for ch in row])
    return Grid(cells, focus)


if __name__ == "__main__":
    # ── 1. Classic Game of Life: Blinker oscillator ─────
    print("=" * 50)
    print("GAME OF LIFE — BLINKER (period-2 oscillator)")
    print("=" * 50)
    print()

    blinker = make_grid([
        ".....",
        ".....",
        ".###.",
        ".....",
        ".....",
    ], focus=(2, 2))

    run(blinker, 4)

    # ── 2. Glider ───────────────────────────────────────
    print("=" * 50)
    print("GAME OF LIFE — GLIDER")
    print("=" * 50)
    print()

    glider = make_grid([
        ".#........",
        "..#.......",
        "###.......",
        "..........",
        "..........",
        "..........",
        "..........",
        "..........",
        "..........",
        "..........",
    ], focus=(0, 1))

    run(glider, 8)

    # ── 3. Co-Kleisli composition demo ──────────────────
    print("=" * 50)
    print("CO-KLEISLI COMPOSITION: life_rule <=< life_rule")
    print("=" * 50)
    print()
    print("Two generations in one composed arrow.")
    print()
    print("To know my state at gen 2, I need my neighbors' states")
    print("at gen 1 — which each required THEIR neighbors at gen 0.")
    print("That's a two-hop computation: co-Kleisli runs life_rule")
    print("everywhere first (building gen 1), then life_rule reads")
    print("from that new grid to produce gen 2.")
    print()

    blinker2 = make_grid([
        ".......",
        ".......",
        ".......",
        ".###...",
        ".......",
        ".......",
        ".......",
    ], focus=(3, 2))

    print("Generation 0:")
    print(blinker2)
    print()

    # One step: extend life_rule
    gen1 = blinker2.extend(life_rule)
    print("Generation 1 (extend life_rule):")
    print(gen1)
    print()

    # Two steps the obvious way: extend twice
    gen2 = gen1.extend(life_rule)
    print("Generation 2 (extend life_rule again):")
    print(gen2)
    print()

    # Same thing via co-Kleisli composition: life_rule <=< life_rule
    two_step = co_kleisli_compose(life_rule, life_rule)
    gen2_composed = blinker2.extend(two_step)
    print("Generation 2 (co-Kleisli: life_rule <=< life_rule):")
    print(gen2_composed)
    print()

    print("Same result — but the composed arrow is a single")
    print("Grid[Cell] -> Cell function that jumps two generations.")
