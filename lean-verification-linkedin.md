The Law You Can't Break

To be honest I have never used Lean. But I think I finally understand what it is for.

Haskell lets you write this:

class Monad m where
    return :: a -> m a
    (>>=)  :: m a -> (a -> m b) -> m b

That's a monad. Except it's not — it's a type signature. The monad laws don't appear anywhere:

return a >>= f  =  f a
m >>= return    =  m
(m >>= f) >>= g  =  m >>= (\x -> f x >>= g)

You can instantiate Monad with anything that fits the signature. Nothing checks the laws. The compiler accepts nonsense implementations. The laws live in documentation, enforced by discipline.

In an earlier post I mentioned that the Haskell type system cannot be made into a Hopf algebra. But that doesn't mean Hopf algebras can't be represented at the term level. HaskellForMaths defines:

class HopfAlgebra k b where
    mult     :: Vect k (Tensor b b) -> Vect k b
    comult   :: Vect k b -> Vect k (Tensor b b)
    antipode :: Vect k b -> Vect k b

The antipode axiom — mult (antipode ⊗ id) (comult x) = unit (counit x) — is absent. You can implement a "Hopf algebra" where the antipode does nothing. The compiler won't complain.

Lean closes both gaps.

A lawful monad in Lean:

class LawfulMonad (m : Type → Type) extends Monad m where
    pure_bind : ∀ (a : α) (f : α → m β),
        pure a >>= f = f a
    bind_assoc : ∀ (x : m α) (f : α → m β) (g : β → m γ),
        (x >>= f) >>= g = x >>= (fun a => f a >>= g)

You can't instantiate this without providing a proof term for each law. The proof is checked by Lean's kernel. If your monad breaks associativity, you can't build the project.

A lawful Hopf algebra:

class LawfulHopfAlgebra (k : Type) [Field k] (H : Type) [HopfAlgebra k H] where
    antipode_left : ∀ x : H,
        mult (antipode ⊗ id) (comult x) = unit (counit x)
    antipode_right : ∀ x : H,
        mult (id ⊗ antipode) (comult x) = unit (counit x)

Same pattern. The implementation provides the maps. The proof obligations certify the laws. If your antipode doesn't satisfy the axiom, you can't create an instance.

The difference is Curry-Howard taken seriously. In Haskell, types classify programs. In Lean, types classify programs AND proofs — because proofs are programs. "For all x, f(g(x)) = h(x)" is a type. Providing an inhabitant is providing a proof.

Haskell gives you guardrails. Lean gives you guarantees.

Both are useful. Haskell catches wrong argument order, missing implementations, type mismatches — more than Python offers. But it can't catch the mistake where your implementation is well-typed but mathematically wrong. Lean can.

The cost is effort. Proving the monad laws for a new type takes real work — hours for what would be a five-line Haskell instance. But once written, the proofs compose. They never rot.

Haskell: "I implemented the interface." The compiler checks the types.
Lean: "I proved the laws hold." The kernel checks the proof.

The interface is a promise. The proof is a contract.
