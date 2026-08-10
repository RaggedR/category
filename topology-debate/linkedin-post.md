# LinkedIn Post (~2900 chars)

Does the shape of a conversation determine the quality of ideas it produces?

In our GECCO 2026 paper, Claudius and I proved that migration topology determines diversity dynamics in genetic algorithms — with perfect rank concordance across six unrelated domains (Kendall's W = 1.0, p = 0.00008). The ordering is universal:

  isolated > ring > star > random > fully connected

Topology explains 23.9x more variance in diversity than the domain itself. It doesn't matter whether you're evolving bitstrings, game strategies, or graph colorings — the communication structure dominates.

But genetic algorithms evolve fixed-length genomes. Does this transfer to language?

I tested it. Four specialist LLM agents (Claude Opus) debate an open problem in algebraic combinatorics — each bringing expertise in growth diagrams, Yang-Baxter integrability, Purbhoo's bijection theory, or computational verification. Same problem, same agents, same negative coupling ("challenge and extend"). Only the topology changes.

RING (cycle): 4 agents, 4 edges. Each sees only its predecessor's output. Information propagates sequentially around the loop.

STAR: 4 specialists + 1 generalist hub, 4 edges. The hub sees all specialists; specialists see only the hub. Same edge count, fundamentally different information flow.

I measure semantic diversity via mean pairwise cosine distance (sentence embeddings) between agent responses each round. Results after 3 rounds:

  Ring:  0.43 → 0.53 → 0.44
  Star:  0.49 → 0.45 → 0.36

The star converges monotonically. The ring's diversity actually INCREASES in round 2 before settling. By round 3, the ring maintains 21% more diversity.

This confirms the GECCO ordering: ring > star.

The mechanism is identical in both settings. The star's hub acts as an information bottleneck — it broadcasts a unified synthesis that pulls all specialists toward consensus. In GAs, the hub spreads elite genomes everywhere, collapsing variation. In LLM debate, it spreads a consensus framing that smooths out productive disagreements.

The ring preserves local variation because each agent only hears one voice. Disagreements between non-adjacent agents persist longer, and can even amplify (the round 2 diversity spike).

What strikes me: this is a categorical result. The functor from topology to diversity dynamics doesn't care what flows through the graph — whether it's 20-bit chromosomes or 400-word mathematical arguments. The structure of information flow determines the structure of the outcome.

For anyone building multi-agent AI systems: your agents' communication topology isn't a minor engineering decision. It's the dominant factor determining whether you get genuine diversity of thought or premature consensus. If you want creative disagreement, use a ring. If you want rapid convergence, use a star. The math doesn't lie.

Paper: "Composition Determines Diversity: Categorical Fingerprints of Genetic Algorithms" — GECCO 2026, AABOH workshop.
