"""Partition refinement on the one-loop graph.

One vertex u, one loop a.  A labelling is a sequence l_0, l_1, ...
with l_k = l(a^k).  The table at a^k is the tail (l_k, l_{k+1}, ...).
Classes after n refinement steps = distinct windows of length n+1.
A split at step n = a new class appearing at window length n+1.
"""

import random

L = 100_000  # sequence length; windows use positions 0..L-n-1
DEPTH = 25   # refinement steps


def fibonacci_word(length):
    s = "0"
    while len(s) < length:
        s = s.replace("0", "0x").replace("1", "0").replace("x", "1")
    return s[:length]


def periodic(length, period="0110"):
    return (period * (length // len(period) + 1))[:length]


def random_word(length, seed=42):
    rng = random.Random(seed)
    return "".join(rng.choice("01") for _ in range(length))


def classes(word, n):
    """Number of refinement classes at step n = distinct windows of length n+1."""
    w = n + 1
    return len({word[i:i + w] for i in range(len(word) - w + 1)})


labellings = {
    "periodic 0110": periodic(L),
    "Fibonacci (Sturmian)": fibonacci_word(L),
    "random": random_word(L),
}

header = f"{'step n':>6} | " + " | ".join(f"{name:>20}" for name in labellings)
print(header)
print("-" * len(header))
prev = {name: 1 for name in labellings}
for n in range(DEPTH + 1):
    row = []
    for name, word in labellings.items():
        c = classes(word, n)
        splits = c - prev[name]
        prev[name] = c
        row.append(f"{c:>13} (+{splits})")
    print(f"{n:>6} | " + " | ".join(f"{cell:>20}" for cell in row))
