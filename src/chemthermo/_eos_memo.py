"""A bounded, call-local memo for repeated equation-of-state solves (ADR-0030).

Internal (ADR-0001): nothing here is re-exported from ``chemthermo``.

What it is for
--------------
One ``flash_tp`` or ``stability_tp`` call asks its equation of state for the
*same* density/compressibility roots several times over. ADR-0023 removed the
part of that a single caller could see - the minimum-Gibbs rule asking for two
branches at one composition - with ``ln_fugacity_branches``. What it could not
see is the duplication *between* callers, because no caller knows what another
one already solved. Measured at ``d8814de`` on the ADR-0023 workload, counting
``PCSAFTEOS._density_roots`` calls at exactly equal ``(names, T, P, x)``:

===============================  =======  ==========  ======
case                             solves   distinct    repeat
===============================  =======  ==========  ======
``pcsaft-vle-methane-hexane``        137         102    25.5 %
``pcsaft-lle-water-hexane``          275         238    13.5 %
``pcsaft-polymer-lle``               601         512    14.8 %
===============================  =======  ==========  ======

and the repeats are the three pairs ADR-0023's roadmap named: a
``fugacity_coefficients`` followed by ``phase_identity`` at the same state, a
``fugacity_coefficients`` followed by ``ln_fugacity_branches``, and the
post-split stability test re-evaluating a composition the split just left.

Why the memo is in the call and not on the model
------------------------------------------------
ADR-0023 rejected "cache the last root solve on the model instance", and the
reason still holds: both model classes are frozen dataclasses, a hidden mutable
cache makes them effectively stateful and thread-unsafe, and a stale key in one
would be a *wrong thermodynamic answer* rather than a slow one. So nothing is
written to a model here. The memo lives in the **call**: it is installed by
:func:`scoped` (which :func:`chemthermo.flash_tp` and
:func:`chemthermo.stability_tp` wear as a decorator), it is read through a
:class:`~contextvars.ContextVar`, and it is discarded when the call returns.
A ``ContextVar`` is per-thread and per-async-task, so two concurrent flashes
never see each other's memo, and a model instance shared between them still
carries no state at all.

Why a hit is bit-identical
--------------------------
A hit returns **the identical object** the first solve produced, not a
recomputation of it, so every double downstream of it is the same double. The
key is the solver's own arguments compared for exact equality, which is why
:func:`composition_key` takes the trouble to separate ``-0.0`` from ``0.0``:
those two compare equal and hash equal in Python while being different bit
patterns, and a key that merges them would be the class of stale key ADR-0023
was worried about. Nothing else in a key can compare equal without being the
same double (``nan`` never compares equal to itself, so a ``nan`` argument
simply misses and is recomputed).

The bound is a memory bound, never a numerical one: evicting an entry can only
cause the solver to redo a solve that produces the same doubles again, so no
result anywhere depends on ``max_entries`` or on the eviction order.

Instrumentation
---------------
The memo counts ``hits`` and ``misses``. They are deliberately **not** put in
any result's ``diagnostics``: that mapping is compared whole by the bit-identity
fixture, and a new key there is a changed answer by this repository's own rule.
A caller that wants the counters opens the scope itself - ``scoped`` reuses an
already-active memo rather than nesting a second one, so::

    with chemthermo._eos_memo.activated() as memo:
        chemthermo.flash_tp(...)
    memo.hits, memo.misses

is how ``tests/test_eos_memo.py`` reads them. Note that widening the scope that
way also widens the memo: it then spans every call inside the ``with``, which
is still exact (the key is exact) but is no longer the per-call lifetime the
library itself uses.
"""

from __future__ import annotations

import functools
import math
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Any, Callable, Hashable, Iterator, ParamSpec, Sequence, TypeVar

#: Entries kept before the oldest is evicted. Four thousand keys is far more
#: than any state in this repository reaches (the largest measured is 512, the
#: polymer case) while bounding one call's memo to a few megabytes for a model
#: whose stored value is a handful of floats. See the module docstring: the
#: bound cannot move a number, only the number of solves.
MAX_ENTRIES = 4096

#: Returned by :meth:`EosSolveMemo.lookup` when the key is absent. A sentinel
#: rather than ``None`` because ``None`` is a legitimate stored value.
MISS: Any = object()

_POSITIVE_ZERO = "+0.0"
_NEGATIVE_ZERO = "-0.0"


class EosSolveMemo:
    """A bounded FIFO memo of one call's equation-of-state solves.

    Attributes:
        hits: Lookups that were served from the memo.
        misses: Lookups that were not.
        max_entries: The eviction bound (see :data:`MAX_ENTRIES`).
    """

    __slots__ = ("_entries", "hits", "max_entries", "misses")

    def __init__(self, max_entries: int = MAX_ENTRIES) -> None:
        if max_entries < 1:
            raise ValueError("max_entries must be at least 1.")
        #: key -> (owner, value). The owner is kept alive so that a key built
        #: from ``id(model)`` cannot be aliased by a later object reusing that
        #: id; see :meth:`store`.
        self._entries: dict[Hashable, tuple[object, Any]] = {}
        self.max_entries = max_entries
        self.hits = 0
        self.misses = 0

    def __len__(self) -> int:
        return len(self._entries)

    def lookup(self, key: Hashable) -> Any:
        """The stored value for ``key``, or :data:`MISS`."""
        entry = self._entries.get(key)
        if entry is None:
            self.misses += 1
            return MISS
        self.hits += 1
        return entry[1]

    def store(self, key: Hashable, value: Any, owner: object = None) -> None:
        """Remember ``value`` under ``key``, evicting the oldest entry if full.

        ``owner`` is the object whose ``id`` the key was built from, if any. It
        is stored purely to hold a strong reference: while the entry lives the
        object cannot be collected, so its ``id`` cannot be handed to a
        different object that would then read this entry as its own.
        """
        entries = self._entries
        if key not in entries and len(entries) >= self.max_entries:
            entries.pop(next(iter(entries)))
        entries[key] = (owner, value)


_ACTIVE: ContextVar[EosSolveMemo | None] = ContextVar("chemthermo_eos_solve_memo", default=None)


def active_memo() -> EosSolveMemo | None:
    """The memo of the call in progress, or ``None`` outside one.

    A model consults this: outside a ``flash_tp`` / ``stability_tp`` call there
    is no memo and every solve runs exactly as it did before ADR-0030.
    """
    return _ACTIVE.get()


@contextmanager
def activated(memo: EosSolveMemo | None = None) -> Iterator[EosSolveMemo]:
    """Run the block with a call-local memo installed, and yield it.

    Re-entrant by design: a ``stability_tp`` call made *inside* a ``flash_tp``
    call finds the flash's memo already active and keeps it, so one flash is
    one memo rather than one per nested stability test. The installer is the
    only one that resets the context variable, so the memo's lifetime is
    exactly the outermost scope's.

    ``memo`` supplies the object to install instead of a default-sized one,
    which is how ``tests/test_eos_memo.py`` exercises the eviction bound. It is
    ignored when a memo is already active, for the same re-entrancy reason.
    """
    existing = _ACTIVE.get()
    if existing is not None:
        yield existing
        return
    installed = memo if memo is not None else EosSolveMemo()
    token = _ACTIVE.set(installed)
    try:
        yield installed
    finally:
        _ACTIVE.reset(token)


_P = ParamSpec("_P")
_R = TypeVar("_R")


def scoped(function: Callable[_P, _R]) -> Callable[_P, _R]:
    """Decorator form of :func:`activated`, for the two public entry points."""

    @functools.wraps(function)
    def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        with activated():
            return function(*args, **kwargs)

    return wrapper


def composition_key(composition: Sequence[float]) -> tuple[object, ...]:
    """A hashable key for ``composition`` that never merges two different bit patterns.

    ``tuple(composition)`` would be the obvious key and is wrong in exactly one
    way: ``-0.0 == 0.0`` and ``hash(-0.0) == hash(0.0)``, so a composition
    carrying a negative zero would read an entry solved for a positive one.
    Nothing measured here produces a negative zero and both models would almost
    certainly return the same doubles for either - but "almost certainly" is
    not the standard a memo whose miss is a wrong answer is held to, so the two
    zeros are given different keys and the question does not arise.
    """
    return tuple(
        value if value else (_POSITIVE_ZERO if math.copysign(1.0, value) > 0.0 else _NEGATIVE_ZERO)
        for value in composition
    )
