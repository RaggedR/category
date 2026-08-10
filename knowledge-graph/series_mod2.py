"""Refinement class counts for coefficient sequences mod 2.

One-loop graph, values in F_2: label a^k with the k-th coefficient
mod 2.  Classes at step n = distinct windows of length n+1.

  Fibonacci numbers mod 2 : rational series      -> expect bounded
  Catalan numbers mod 2   : algebraic (Christol) -> expect linear
  Partition numbers mod 2 : (conjecturally) transcendental
"""

L = 200_000
DEPTH = 22


def fib_mod2(length):
    a, b, out = 0, 1, []
    for _ in range(length):
        out.append(a & 1)
        a, b = b, a + b  # keep mod 2 to stay fast
        a %= 2
        b %= 2
    return "".join(map(str, out))


def catalan_mod2(length):
    # C_n odd  iff  n = 2^k - 1
    return "".join("1" if (n + 1) & n == 0 else "0" for n in range(length))


def partition_mod2(length):
    # Euler's pentagonal recurrence, all arithmetic mod 2
    p = [0] * length
    p[0] = 1
    for n in range(1, length):
        total, k = 0, 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > n:
                break
            total += p[n - g1]
            if g2 <= n:
                total += p[n - g2]
            k += 1
        p[n] = total & 1
    return "".join(map(str, p))


words = {
    "Fibonacci mod 2 (rational)": fib_mod2(L),
    "Catalan mod 2 (algebraic)": catalan_mod2(L),
    "partition mod 2 (transcendental?)": partition_mod2(L),
}

for name, w in words.items():
    print(f"\n=== {name}   starts: {w[:30]}...")
    prev = 1
    for n in range(DEPTH + 1):
        size = n + 1
        c = len({w[i:i + size] for i in range(len(w) - size + 1)})
        print(f"  n={n:>2}  classes={c:>7}  (+{c - prev})")
        prev = c
