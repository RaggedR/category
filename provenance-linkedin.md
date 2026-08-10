Don't Return a Number. Return the Proof.

I've been studying an AI-assisted tax platform where the tax engine doesn't output a final liability. It outputs a directed acyclic graph. And the mechanism that makes this work is the same one that scares people away from Haskell: a monad.

Every node in the graph records which rule fired, which inputs it consumed, and which value it produced. The final tax figure is just the root node. The SA100 form filed with HMRC, the audit trail, and the human-readable explanation ("why is box 1 £42,500?") are all views of this same graph. They can't disagree -- they're derived from the same data structure.

The mechanism that threads the graph through the computation is a Writer monad. As each tax rule fires, it appends a node to the trace. The programmer never passes the trace manually -- the monad carries it. Without it, every function would take a trace and return a trace. You could forget. The monad makes that impossible.

Now add AI. The system uses LLMs for receipt OCR and expense classification. Models hallucinate. A receipt says £420 but the LLM extracts £4,200. No type system catches this -- it's a valid amount.

So how do you catch it? Every value carries provenance as a field on the value itself -- not in a log table, not in a sidecar database. A Fact like "office supplies: £4,200" carries a reference to the source document, byte offsets, extraction method, and confidence score. This threads through every transformation via the same monadic machinery.

The return says £4,200. The user clicks through the provenance link and sees the original receipt says £420. The error wasn't prevented. But it cannot hide.

The document puts it simply: "A recommendation made on Tuesday must survive a challenge on Friday -- by the user, by HMRC, or by a compliance reviewer." That's what defensible means. Not that every number is right, but that every number can be interrogated. The system doesn't say "the code computed it." It produces a chain of rules, inputs, and source documents an auditor can verify independently.

The monad from a Haskell "hello world" tutorial is the same mechanism threading an audit trail through every tax computation. That's not accidental complexity. That's the product.

Don't return a number. Return the proof.
