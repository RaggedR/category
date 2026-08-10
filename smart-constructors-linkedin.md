The Bouncer at the Door

I've been studying an AI-assisted tax platform with an unusual architecture: a Python outer shell wrapping a strongly-typed Haskell inner core. The most important design decision isn't about either language. It's about what happens at the boundary between them.

The system uses AI for receipt OCR and expense classification. Models hallucinate. A scanned receipt might say £420 but the LLM extracts £4,200. A bank statement credit might get classified as self-employment income when it's a refund. These aren't edge cases -- they're Tuesday.

So the architects made a rule: AI outputs never bypass the core. Every value the AI produces -- a classified expense, an extracted figure, a generated explanation -- must pass through a "smart constructor" before it can influence any computation.

A smart constructor is a function that is the only way to create a value of a given type. The plain data constructor is hidden. You can't just wrap a raw number and call it a TaxYear. You call mkTaxYear, which checks that the year is valid and returns either a verified TaxYear or an error. If it returns a TaxYear, that value is trustworthy for the rest of its life. No downstream code needs to re-check it, because there is no other way it could have been created.

In Haskell, this is enforced by the module system. The module exports the type name (so other code can use it in signatures) but not the data constructor (so other code can't fabricate values). In Python, you could write a factory method and hide the constructor by convention -- but nothing stops someone bypassing it. Convention breaks under deadline pressure. Module boundaries don't.

The tax platform applies this at every point where data crosses from the untrusted outer shell into the trusted core. A raw extraction from an LLM is just a string with a confidence score. It becomes a Fact -- a typed, validated, provenance-carrying value -- only by passing through a smart constructor that checks the amount, the category, the tax year, and attaches a reference back to the source document. If the LLM hallucinated, the constructor rejects it. The hallucination never reaches the tax engine.

This creates a powerful invariant: if a value exists inside the core, it is valid. Not "it was probably validated." Not "it should have been checked." It is valid, because there is no other way it could have come into existence.

The practical payoff is that every function inside the core can trust its inputs. No defensive re-checking. No "just in case" validation scattered through the codebase. The smart constructor at the boundary did the hard work once, and the type carries the proof everywhere the value goes.

The document calls this "correct by construction." A TaxYear tagged as 2025 is correct because the only function that can produce one verified that it was 2025. You don't validate before use because it could not have been created wrong.

This matters beyond tax. Every system that takes input from an unreliable source -- an LLM, a user form, an external API, a CSV upload -- faces the same question: where do you draw the line between "untrusted data" and "trusted data"? The smart constructor pattern says: draw it once, draw it sharp, and make it impossible to cross without validation. Then stop worrying.

The best code isn't code that handles every possible error. It's code that structurally cannot receive invalid input in the first place.
