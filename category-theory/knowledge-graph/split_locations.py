"""Where do the splits happen?

For a rotation word with slope alpha, the class of a window w (length n)
splits at step n when both w0 and w1 occur ("right-special" window).
Sturmian: exactly one right-special window per length.  We track that
window and how it grows, for quadratic and non-quadratic slopes.
"""

import math

L = 200_000
DEPTH = 20

ALPHAS = {
    "Fibonacci  a=1/phi^2": (math.sqrt(5) - 1) / 2 ** 2 if False else 1 / ((1 + math.sqrt(5)) / 2) ** 2,
    "Pell       a=sqrt2-1": math.sqrt(2) - 1,
    "non-quadr  a=1/e":     1 / math.e,
}


def rotation_word(alpha, length):
    return "".join(
        str(math.floor((k + 1) * alpha) - math.floor(k * alpha))
        for k in range(length)
    )


for name, alpha in ALPHAS.items():
    word = rotation_word(alpha, L)
    print(f"\n=== {name}   word: {word[:40]}...")
    specials = []
    for n in range(1, DEPTH + 1):
        seen = {}
        for i in range(len(word) - n):
            w = word[i:i + n]
            seen.setdefault(w, set()).add(word[i + n])
        sp = [w for w, nxt in seen.items() if len(nxt) == 2]
        specials.append(sp)
        shown = sp[0] if len(sp) == 1 else sp
        print(f"  n={n:>2}  #special={len(sp)}  split window: {shown}")

    # If each length has exactly one special window, check the chain
    if all(len(sp) == 1 for sp in specials):
        chain_ok = all(
            specials[n][0].endswith(specials[n - 1][0])
            for n in range(1, len(specials))
        )
        # prepended letters: w_{n+1} = c_{n+1} . w_n ; c_1 = w_1
        c = specials[0][0] + "".join(
            specials[n][0][0] for n in range(1, len(specials))
        )
        print(f"  suffix chain: {chain_ok}")
        print(f"  prepended letters c : {c}")
        print(f"  word prefix         : {word[:len(c)]}")
        print(f"  c == word prefix? {c == word[:len(c)]}")
        rev = specials[-1][0][::-1]
        print(f"  reversed special == word prefix? {rev == word[:len(rev)]}")
