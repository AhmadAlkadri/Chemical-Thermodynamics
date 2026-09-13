"""Shared binary-interaction-parameter (``kij``) plumbing.

Internal module (ADR-0001: anything outside the public API sources is
internal). It holds the canonicalization and matrix-building helpers that
``PengRobinsonEOS`` introduced in ADR-0006, so that a second model family
(PC-SAFT, ADR-0014) can reuse the *same* name-keyed, order-agnostic contract
instead of growing a parallel one.

Nothing here changed when the code moved out of
``chemthermo/models/peng_robinson.py``: the functions are byte-for-byte the
former ``_canonicalize_kij`` and ``PengRobinsonEOS._kij_matrix`` bodies, with
the latter taking a sequence of component names instead of a ``Mixture`` so it
can serve a model that is not built around ``Mixture``. Peng-Robinson results
are therefore bit-identical (``tests/test_pr_eos.py``,
``tests/test_flash_refactor_bit_identity.py``).
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from ..data import normalize_name
from ..exceptions import ModelError

#: Canonical stored form of a per-pair kij matrix: a sorted tuple of
#: ``((name_a, name_b), value)`` entries with ``name_a < name_b`` (both
#: normalized via ``chemthermo.data.normalize_name``). Kept immutable so a
#: frozen dataclass stays hashable/comparable and its repr is deterministic.
KijPairs = tuple[tuple[tuple[str, str], float], ...]
KijInput = float | Mapping[tuple[str, str], float]


def canonicalize_kij(kij: KijInput | KijPairs, *, model: str) -> float | KijPairs:
    """Normalize constructor input for a model's ``kij`` field.

    A scalar is returned unchanged (as a ``float``); it is applied to every
    off-diagonal pair and never to the diagonal. A mapping from unordered
    component-name pairs to values is normalized (name case/whitespace via
    ``normalize_name``, pair order) into a sorted tuple. Both orders of a pair
    must agree if both are given; a pair naming the same component twice is
    rejected. Values are looked up per-mixture later, so pair names unknown to
    any particular mixture are simply never used.

    Idempotent: re-running this on an already-canonical tuple (as happens on
    ``dataclasses.replace``) returns it unchanged rather than re-validating,
    since it was already validated when first constructed.

    ``model`` only names the caller in error messages.
    """
    if isinstance(kij, tuple):
        return kij
    if isinstance(kij, Mapping):
        canonical: dict[tuple[str, str], float] = {}
        for (name_a, name_b), value in kij.items():
            key_a = normalize_name(name_a)
            key_b = normalize_name(name_b)
            if not key_a or not key_b:
                raise ModelError(f"{model} kij component names must be non-empty.")
            if key_a == key_b:
                raise ModelError(
                    f"{model} kij pair components must be distinct "
                    f"(got {name_a!r} paired with itself)."
                )
            pair_key = (key_a, key_b) if key_a < key_b else (key_b, key_a)
            value_f = float(value)
            if pair_key in canonical and canonical[pair_key] != value_f:
                raise ModelError(
                    f"Conflicting kij values given for pair {pair_key!r}: "
                    f"{canonical[pair_key]!r} vs {value_f!r}."
                )
            canonical[pair_key] = value_f
        return tuple(sorted(canonical.items()))
    return float(kij)


def kij_matrix(kij: KijInput | KijPairs, names: Sequence[str]) -> np.ndarray:
    """Build the dense ``n x n`` kij matrix for ``names``' order.

    The diagonal is always zero. A scalar ``kij`` fills every off-diagonal
    entry; a per-pair mapping fills only the pairs it names (by normalized
    component name) and defaults missing pairs to zero.
    """
    n = len(names)
    matrix = np.zeros((n, n), dtype=float)

    if isinstance(kij, (int, float)):
        if kij != 0.0:
            matrix[:, :] = kij
            np.fill_diagonal(matrix, 0.0)
        return matrix

    # ``kij`` is always the canonical KijPairs tuple here (never a raw Mapping
    # at runtime; the owning dataclass's __post_init__ already converted it).
    # Both branches call the same dict(...) constructor -- the isinstance split
    # exists only so pyright resolves a single dict() overload per branch
    # instead of over-widening a `Mapping | KijPairs` union.
    pairs: dict[tuple[str, str], float]
    if isinstance(kij, Mapping):
        pairs = dict(kij)
    else:
        pairs = dict(kij)
    canonical_names = [normalize_name(name) for name in names]
    for i in range(n):
        for j in range(i + 1, n):
            key = (
                (canonical_names[i], canonical_names[j])
                if canonical_names[i] < canonical_names[j]
                else (canonical_names[j], canonical_names[i])
            )
            value = pairs.get(key, 0.0)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix
