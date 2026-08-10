# Trees of Trees

## S-Expressions

An **S-expression** over a set of atoms A is a labelled rose tree: either an atom a ∈ A (a leaf), or a node labelled a ∈ A with a list of S-expression children.

SExpr(A) = μT. A × List(T)

---

## The Problem

A **Lisp interpreter written in Lisp** is an S-expression that takes another S-expression as input and evaluates it. Both the interpreter and the program it interprets are trees of the same type. The interpreter receives the program as quoted data — a tree embedded inside a tree.

### Part 1: Trees of trees

The Lisp expression `(eval '(+ 1 2))` consists of an interpreter's `eval` call (the outer tree) containing a quoted program `(+ 1 2)` (the inner tree). Decompose it as an element of T(T(X)) — a tree whose leaves are trees.

### Part 2: Writing a Lisp interpreter in Lisp is self-composition

A Lisp interpreter written in Lisp is an element I ∈ T. The program it interprets is another element P ∈ T. The combined expression (the interpreter applied to the quoted program) is an element of T(T(X)). Explain why this is **container self-composition** — the same container on both sides — and why this requires homoiconicity (code and data share the same syntax).

A Lisp interpreter interpreting a Lisp interpreter interpreting a program is an element of T(T(T(X))) = T³(X).

### Part 3: Folding the tree of trees

A **catamorphism** (fold) against an algebra α : A × List(V) → V processes a tree bottom-up: each node sees only the results of its children, never their structure. The tree is consumed and a value is produced.

The Lisp evaluator's algebra is:

```
α(n, [])           = n                   (integer literal)
α(+, [v₁, v₂])    = v₁ + v₂             (addition)
α(×, [v₁, v₂])    = v₁ × v₂             (multiplication)
α(if, [v₁, v₂, v₃]) = v₂ if v₁ ≠ 0, else v₃
α(λ, [body])       = a closure
α(f, [v₁, …, vₙ])  = apply f to arguments
```

The algebra pattern-matches on labels. This is where the semantics lives — `+` means addition because α says so.

The monad multiplication μ : T(T(X)) → T(X) is a different operation: it **grafts** inner trees into outer leaves, producing a single tree. It does not look at labels. It rearranges; it does not interpret.

Show that eval is a catamorphism, not monad multiplication.

### Part 4: The theorem

**Theorem.** *Writing a Lisp interpreter in Lisp decomposes into three components:*

1. ***Self-composition of labelled rose trees.*** *The interpreter applied to a program is an element of T ∘ T — the syntax container composed with itself. The outer tree is the interpreter; the inner trees are quoted program fragments.*

2. ***Quote is the constructor.*** *Quote places non-trivial inner trees at leaves of the outer tree. Without quote, all inner trees are trivial and the composition collapses to a single tree.*

3. ***Eval is a fold.*** *The interpreter evaluates the program by folding T(T(X)) bottom-up against an algebra α. This is a catamorphism — it consumes the tree and produces values. It is not monad multiplication, which produces trees.*
