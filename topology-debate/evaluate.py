"""Evaluation: diversity metrics and cross-domain connection detection."""

from itertools import combinations

import numpy as np
from sentence_transformers import SentenceTransformer

from state import RING_ORDER


class DiversityMeter:
    """Measures semantic diversity via pairwise cosine distance.

    Reuses the pattern from ~/git/orchestration/experiment.py.
    Embeddings are L2-normalized, so cosine distance = 1 - dot product.
    """

    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        self.model = SentenceTransformer(model_name)

    def measure(self, texts: list[str], labels: list[str] | None = None) -> dict:
        """Compute pairwise diversity metrics for a set of texts."""
        if len(texts) < 2:
            return {"mean_distance": 0.0, "per_pair": {}}

        emb = self.model.encode(texts, normalize_embeddings=True)
        sim = emb @ emb.T
        dist = 1.0 - sim

        pairs = list(combinations(range(len(texts)), 2))
        distances = [float(dist[i, j]) for i, j in pairs]

        if labels is None:
            labels = [str(i) for i in range(len(texts))]

        per_pair = {
            f"{labels[i]}<->{labels[j]}": d
            for (i, j), d in zip(pairs, distances)
        }

        return {
            "mean_distance": float(np.mean(distances)),
            "min_distance": float(np.min(distances)),
            "max_distance": float(np.max(distances)),
            "per_pair": per_pair,
        }


# -- Cross-domain connection detection ----------------------------------------

DOMAIN_KEYWORDS = {
    "growth_diagram": {
        "fomin", "burge", "face label", "local rule", "recording",
        "growth diagram", "cylindric", "langer",
    },
    "yang_baxter": {
        "yang-baxter", "r-matrix", "commuting", "free-fermion", "free fermion",
        "vertex weight", "spectral parameter", "integrability", "transfer matrix",
    },
    "bijection": {
        "purbhoo", "mosaic", "migration", "jeu de taquin", "jdt",
        "schutzenberger", "reading order",
    },
    "computational": {
        "enumerate", "lr_puzzles", "branching", "verify", "counterexample",
        "c^ν", "minimal example", "code",
    },
}


def detect_cross_domain_connections(responses_in_round: list[dict]) -> list[str]:
    """Find terms appearing in agent X's output that come from agent Y's domain.

    This is a heuristic proxy for "novel cross-domain insight" — when an agent
    uses another domain's vocabulary, it suggests conceptual bridging.
    """
    connections = []
    for entry in responses_in_round:
        agent_id = entry["agent_id"]
        text_lower = entry["text"].lower()
        for other_id, keywords in DOMAIN_KEYWORDS.items():
            if other_id == agent_id:
                continue
            found = [kw for kw in keywords if kw in text_lower]
            if len(found) >= 2:
                connections.append(
                    f"{entry['agent_name']} borrowed from {other_id}: {found}"
                )
    return connections
