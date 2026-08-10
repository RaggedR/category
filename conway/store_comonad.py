"""
Co-Kleisli composition for parameterized tests using the Store comonad.

The Store comonad models a collection of values indexed by position,
with a focused position.  A co-Kleisli arrow  Store s a → b  is a
test that reads from the ambient environment (including peeking at
other positions) and produces a result.

Composing two such tests means: run the first across all positions
(extend), building a new landscape of results, then the second test
reads from that landscape — comparing results across positions.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, TypeVar

A = TypeVar("A")
B = TypeVar("B")
S = TypeVar("S")


@dataclass
class Store:
    """Store s a = (s → a, s) — values indexed by position, with a focus."""

    lookup: dict          # s → a  (using a dict for finite stores)
    focus: str            # the current position

    def extract(self) -> any:
        """Return the value at the focused position."""
        return self.lookup[self.focus]

    def duplicate(self) -> Store:
        """
        Store s a → Store s (Store s a)

        At every position, place the whole store focused there.
        """
        stores = {}
        for key in self.lookup:
            stores[key] = Store(self.lookup, key)
        return Store(stores, self.focus)

    def extend(self, f: Callable[[Store], any]) -> Store:
        """
        (Store s a → b) → Store s a → Store s b

        Apply f at every position, each time focused there.
        """
        results = {}
        for key in self.lookup:
            refocused = Store(self.lookup, key)
            results[key] = f(refocused)
        return Store(results, self.focus)

    def peek(self, key: str) -> any:
        """Read a value at a different position."""
        return self.lookup[key]

    def keys(self):
        return self.lookup.keys()

    def __repr__(self) -> str:
        lines = []
        for key in sorted(self.lookup.keys()):
            marker = " ←" if key == self.focus else ""
            lines.append(f"  {key}: {self.lookup[key]}{marker}")
        return "Store(\n" + "\n".join(lines) + "\n)"


def co_kleisli_compose(f, g):
    """g <=< f = g ∘ extend(f)"""
    def composed(store: Store) -> any:
        return g(store.extend(f))
    return composed


# ── Example: parameterized database tests ───────────────

def make_test_env() -> Store:
    """
    A test environment: three database configs, each with
    a table of users.  Focused on 'postgres'.
    """
    return Store(
        lookup={
            "sqlite": {
                "users": [
                    {"name": "alice", "email": "alice@test.com", "age": 30},
                    {"name": "bob", "email": "bob@test.com", "age": 25},
                ],
                "config": {"case_sensitive": False, "max_connections": 1},
            },
            "postgres": {
                "users": [
                    {"name": "Alice", "email": "alice@test.com", "age": 30},
                    {"name": "Bob", "email": "bob@test.com", "age": 25},
                ],
                "config": {"case_sensitive": True, "max_connections": 100},
            },
            "mysql": {
                "users": [
                    {"name": "alice", "email": "alice@test.com", "age": 30},
                    {"name": "bob", "email": "bob@test.com", "age": 25},
                ],
                "config": {"case_sensitive": False, "max_connections": 50},
            },
        },
        focus="postgres",
    )


# ── Co-Kleisli arrows (tests) ──────────────────────────

def find_user_by_name(env: Store) -> str | None:
    """
    Store Config DB → str | None

    Search for 'alice' (lowercase) in the focused database.
    Returns the email if found, None otherwise.
    This is a co-Kleisli arrow: it reads from the ambient
    environment (one database config) and returns a single result.
    """
    db = env.extract()
    for user in db["users"]:
        if user["name"].lower() == "alice":
            return user["email"]
    return None


def results_agree(env: Store) -> bool:
    """
    Store Config Result → bool

    Check whether ALL configs produced the same result as
    the focused config.  This PEEKS at other positions —
    it needs the ambient landscape of results, not just its own.

    This is where co-Kleisli composition earns its keep.
    Plain function composition would only give us one result.
    We need extend to run the first test everywhere, building
    the landscape that this test reads from.
    """
    my_result = env.extract()
    for key in env.keys():
        if env.peek(key) != my_result:
            return False
    return True


def find_discrepancies(env: Store) -> list[tuple[str, str]]:
    """
    Store Config Result → list of (config, result) that disagree
    with the focused config.

    Peeks at every config's result to find inconsistencies.
    """
    my_result = env.extract()
    mismatches = []
    for key in env.keys():
        other = env.peek(key)
        if other != my_result:
            mismatches.append((key, other))
    return mismatches


# ── Demo ────────────────────────────────────────────────

if __name__ == "__main__":
    env = make_test_env()

    print("=" * 55)
    print("TEST ENVIRONMENT (Store Config DB, focused on postgres)")
    print("=" * 55)
    print()
    for config in sorted(env.lookup.keys()):
        marker = " ← focused" if config == env.focus else ""
        db = env.lookup[config]
        names = [u["name"] for u in db["users"]]
        print(f"  {config}: users = {names}{marker}")
    print()

    # ── Step 1: extend(find_user_by_name) ───────────────
    print("=" * 55)
    print("STEP 1: extend(find_user_by_name)")
    print("Run the test at every config, building a Store of results")
    print("=" * 55)
    print()

    results = env.extend(find_user_by_name)
    print(results)
    print()
    print("Note: all three configs found alice's email.")
    print("The search was case-insensitive, so even postgres")
    print("(which stores 'Alice') matched.")
    print()

    # ── Step 2: results_agree reads from the landscape ──
    print("=" * 55)
    print("STEP 2: results_agree (reads the landscape)")
    print("Peek at all configs' results — do they agree?")
    print("=" * 55)
    print()

    verdict = results_agree(results)
    print(f"  All configs agree? {verdict}")
    print()

    # ── Composed: results_agree <=< find_user_by_name ───
    print("=" * 55)
    print("COMPOSED: results_agree <=< find_user_by_name")
    print("= results_agree ∘ extend(find_user_by_name)")
    print("=" * 55)
    print()

    composed = co_kleisli_compose(find_user_by_name, results_agree)

    # Run the composed test — it goes straight from DB env to verdict
    print(f"  Composed result: {composed(env)}")
    print()
    print("  Same answer — but now it's a single arrow")
    print("  Store Config DB → bool")
    print()

    # ── Now introduce a bug: case-sensitive search ──────
    print("=" * 55)
    print("NOW WITH A BUG: case-sensitive search")
    print("=" * 55)
    print()

    def find_user_case_sensitive(env: Store) -> str | None:
        """Search for 'alice' (exact match) — will fail on postgres."""
        db = env.extract()
        for user in db["users"]:
            if user["name"] == "alice":  # no .lower()!
                return user["email"]
        return None

    results_buggy = env.extend(find_user_case_sensitive)
    print("Results of case-sensitive search across configs:")
    print(results_buggy)
    print()

    # Compose with find_discrepancies to see what disagrees
    composed_debug = co_kleisli_compose(find_user_case_sensitive, find_discrepancies)
    print(f"  Discrepancies (vs focused config postgres):")
    for config, result in composed_debug(env):
        print(f"    {config} got: {result}")
    print()
    print("  postgres returned None (stores 'Alice', not 'alice')")
    print("  but sqlite and mysql found her — the test is")
    print("  database-dependent. Co-Kleisli composition exposed it.")
