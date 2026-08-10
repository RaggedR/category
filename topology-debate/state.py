"""Typed state and agent configuration for ring and star topology debates."""

import operator
from typing import Annotated, TypedDict


# -- Agent Definitions ---------------------------------------------------------

AGENTS = [
    {
        "id": "growth_diagram",
        "name": "Growth Diagram Agent",
        "system_prompt": (
            "You are a specialist in Fomin growth diagrams, Langer's thesis (Section 4.2), "
            "and the Burge local rule. You understand how face labels encode recording data "
            "in RSK-type correspondences via local rules on directed graphs. Your focus: "
            "whether KTW puzzle branching choices form a growth diagram — the deterministic "
            "vs branching distinction in triangle tiles, and whether the cylindric growth "
            "diagram construction (recursive L_i composition) applies here."
        ),
        "rag_seed_query": "growth diagram local rule Fomin Burge face labels recording tableau",
        "predecessor": "computational",
    },
    {
        "id": "yang_baxter",
        "name": "Yang-Baxter Agent",
        "system_prompt": (
            "You are a specialist in the Yang-Baxter equation, free-fermionic vertex models, "
            "and Zinn-Justin's integrable tilings framework. You understand commuting transfer "
            "matrices T_+ and T_-, R-matrices, and why commutativity guarantees well-defined "
            "growth diagrams. Your focus: the vertex model structure, spectral parameters, "
            "and how free-fermion weights restrict partition differences to <= 1 per row."
        ),
        "rag_seed_query": "Yang-Baxter equation free fermion vertex model transfer matrix commuting",
        "predecessor": "growth_diagram",
    },
    {
        "id": "bijection",
        "name": "Bijection Agent",
        "system_prompt": (
            "You are a specialist in Purbhoo's bijection (puzzle -> mosaic -> migration -> "
            "LR-tableau), jeu de taquin, and Schutzenberger involution. You understand "
            "migration as a sequence of jdt slides and how the reading order of tableau "
            "entries corresponds to an order on puzzle tiles. Your focus: whether growth "
            "diagram face labels can replicate the migration bijection directly, bypassing "
            "the mosaic detour."
        ),
        "rag_seed_query": "Purbhoo mosaic migration jeu de taquin LR tableau bijection puzzle",
        "predecessor": "yang_baxter",
    },
    {
        "id": "computational",
        "name": "Computational Agent",
        "system_prompt": (
            "You are a computational mathematics specialist. You reason about enumerating "
            "puzzle witnesses via lr_puzzles(lam, mu, nu), counting free nabla choices, "
            "and verifying conjectures numerically. Your focus: the minimal non-trivial "
            "example c^(3,2,1)_{(2,1),(2,1)} = 2 with 12 free choices projecting to 3 "
            "tableau entries, and what computational experiments would confirm or refute "
            "the growth diagram conjecture."
        ),
        "rag_seed_query": "puzzle witness transfer matrix branching nabla choice enumerate verify",
        "predecessor": "bijection",
    },
]

VANILLA_AGENT = {
    "id": "vanilla",
    "name": "Vanilla Agent",
    "system_prompt": (
        "You are a mathematician with broad knowledge but no particular specialism. "
        "You synthesize ideas from combinatorics, algebra, and geometry. You are "
        "skeptical of over-specialization and look for the simplest explanation that "
        "unifies multiple viewpoints. Your job is to find common ground, identify "
        "the strongest thread across multiple perspectives, and push for clarity."
    ),
    "rag_seed_query": "Littlewood-Richardson puzzle bijection tableau combinatorics",
}

ALL_AGENTS = AGENTS + [VANILLA_AGENT]

RING_ORDER = [a["id"] for a in AGENTS]

AGENT_BY_ID = {a["id"]: a for a in ALL_AGENTS}

COUPLING_PROMPTS = {
    "negative": (
        "CRITICAL: You must CHALLENGE and EXTEND the previous agent's reasoning. "
        "Find gaps, propose alternative framings, or identify where their domain "
        "knowledge connects non-obviously to yours. Do NOT simply agree."
    ),
    "positive": (
        "BUILD ON the previous agent's reasoning. Take their strongest insight "
        "and develop it further using your specialist knowledge. Show how your "
        "domain provides additional evidence or a more precise formulation."
    ),
}


# -- State Schema --------------------------------------------------------------

class DebateState(TypedDict):
    """Typed state flowing through the debate graph."""
    # Configuration (set once at start)
    problem_statement: str
    max_rounds: int
    coupling: str  # "negative" or "positive"
    use_tools: bool  # whether agents can execute Python, read/write files

    # Debate progression
    current_round: int
    current_agent_idx: int  # index into RING_ORDER

    # Accumulated data (Annotated with operator.add for parallel fan-in)
    responses: Annotated[list[dict], operator.add]  # [{round, agent_id, agent_name, text, ...}]
    metrics: Annotated[list[dict], operator.add]    # [{round, diversity, connections}]
