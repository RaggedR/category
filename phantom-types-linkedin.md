The Type Parameter That Doesn't Exist

I've been studying an AI-assisted tax platform that uses Haskell for its trusted computational core. One design decision keeps standing out: phantom types.

A phantom type is a type parameter that appears in the type signature but never in the value. It exists at compile time and is erased at runtime. Zero cost.

You define money like this:

newtype Money (currency :: Symbol) = Money Int64

The "currency :: Symbol" declares a type parameter that must be a type-level string -- like "GBP" or "USD". It's the phantom: declared on the left side of the definition but never used in the actual data on the right. At runtime, Money "GBP" and Money "USD" are both just 64-bit integers. But at compile time they're different types.

Now you can write a function:

add :: Money c -> Money c -> Money c

The lowercase "c" is a type variable. It can be any currency, but it must be the same in all three positions. If you try to add Money "GBP" to Money "USD", the compiler refuses. Not a runtime error. Not a test failure. The binary never gets built.

The reading rule: lowercase = variable (must match everywhere it appears), uppercase = concrete. Money c means "money in any currency." Money "GBP" means "pounds, specifically."

The tax platform takes this further, parameterising types by jurisdiction AND tax year:

computeReturn :: Rules y -> TaxInput y -> TaxReturn y

Read this as: "computeReturn takes Rules tagged with tax year y, TaxInput tagged with the same y, and returns a TaxReturn tagged with the same y." The compiler enforces that all three agree. You cannot feed 2024 data to 2025 rules -- the y won't unify. You cannot apply UK allowances to US income. Every code path is covered, including the ones nobody thought to test.

Why does this matter for AI?

The platform uses LLMs for receipt OCR and expense classification. Models hallucinate. So every AI output must pass through a smart constructor before it can influence any computation. The phantom types make that trust boundary sharp -- an LLM can't produce a value tagged UK Y2025 unless it passes validation. The type system doesn't care where data came from.

The principle is "correct by construction": once a value enters the system tagged as UK Y2025, that tag was verified at the point of construction. Downstream code doesn't need to check it again -- it could not have been created wrong. You don't write tests for "did someone mix up the tax years?" because the compiler already makes that impossible.

Python handles the IO-heavy work -- OCR, LLM calls, file parsing. Phantom types exist in Rust (PhantomData) and can be simulated in TypeScript (branded types), but are absent from Python -- which is why this system keeps Python at the edges and Haskell at the core, where the arithmetic must be provably correct.

The best validation is the one you never have to write.
