"""PC-SAFT pure-component parameter storage.

Holds the three pure-component parameters every non-associating PC-SAFT
component needs (Gross & Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244):

- ``m``            - number of segments per chain (dimensionless), given
  directly or, for a polymer, as ``segments_per_g`` (mol/g) times ``MW_g_mol``
  (g/mol) - the ``m/M`` convention of ADR-0022,
- ``sigma_A``      - segment diameter, in Angstrom,
- ``epsilon_k_K``  - segment dispersion energy divided by the Boltzmann
  constant, ``epsilon / k_B``, in K.

An **associating** component carries two more (Gross & Sadowski, Ind. Eng.
Chem. Res. 41 (2002) 5510), in an optional :class:`PCSAFTAssociationRecord`
under the record's ``association`` field (ADR-0018):

- ``kappa_ab``         - the dimensionless effective association volume,
- ``epsilon_ab_k_K``   - the association energy divided by ``k_B``, in K,
- ``na`` / ``nb``      - how many sites of each of the two types the molecule
  carries (``1`` / ``1`` is the 2B scheme of the 2002 paper), and
- ``scheme``           - the scheme's conventional name, documentation only:
  ``na`` and ``nb`` are what the model reads.

A record without an ``association`` block is non-associating, and a parameter
set in which no component associates drives exactly the code path of ADR-0014.

The packaged defaults live in ``data/eos/pcsaft.json`` and carry a
``provenance`` block naming their source and how it was verified. Binary
interaction parameters are **not** stored here: ``k_ij`` is a property of the
model instance (``PCSAFTEOS(kij=...)``, ADR-0006), not of a component. Nor are
*cross*-association parameters: they are computed from the pure-component ones
by the combining rules cited in :mod:`chemthermo.eos._pcsaft_association`.

Follows the ``NRTLParameters`` pattern: a packaged JSON payload with a
``schema_version``, a ``model`` label and a list of records, plus a
``from_records`` constructor for user-supplied values. Unknown keys (top level
or inside a record) are ignored, so documentation such as ``provenance`` or a
per-record ``source`` needs no schema bump.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Iterable, Mapping, Sequence

import numpy as np

from ..data import normalize_name
from ..exceptions import ModelError

_PCSAFT_RESOURCE = ("data", "eos", "pcsaft.json")
_SCHEMA_VERSION = 1


class PCSAFTParameterError(ModelError):
    """Raised when PC-SAFT parameters are missing or invalid."""


@dataclass(frozen=True)
class PCSAFTAssociationRecord:
    """Association parameters for one component (ADR-0018).

    Args:
        kappa_ab: Effective association volume ``kappa^{AB}``, dimensionless
            and positive.
        epsilon_ab_k_K: Association energy divided by the Boltzmann constant,
            ``epsilon^{AB}/k_B``, in K. Positive.
        na: Number of type-A ("proton donor") sites per molecule.
        nb: Number of type-B ("proton acceptor") sites per molecule. Only A-B
            bonding is modelled, so both counts must be at least 1.
        scheme: Conventional name of the site scheme (``"2B"``, ``"3B"``,
            ``"4C"``, ...). **Documentation only** - the model reads ``na`` and
            ``nb``. A ``scheme`` naming a count pair that contradicts ``na`` /
            ``nb`` is rejected for the schemes this class knows.
        source: Free-text provenance.

    Only ``na = nb = 1`` (the 2B scheme, which is what Gross & Sadowski (2002)
    use for water and the 1-alkanols) is **validated** in this slice; other
    counts are accepted and implemented by the same general equations but have
    no cross-check behind them. See ADR-0018.
    """

    kappa_ab: float
    epsilon_ab_k_K: float
    na: float = 1.0
    nb: float = 1.0
    scheme: str = "2B"
    source: str = ""

    #: Site counts implied by the scheme names this class knows, so a typo in
    #: ``scheme`` cannot silently disagree with ``na`` / ``nb``.
    _KNOWN_SCHEMES = {"2B": (1.0, 1.0), "3B": (2.0, 1.0), "4C": (2.0, 2.0)}

    def __post_init__(self) -> None:
        for label, value in (
            ("kappa_ab", self.kappa_ab),
            ("epsilon_ab_k_K", self.epsilon_ab_k_K),
        ):
            if not np.isfinite(value) or value <= 0.0:
                raise PCSAFTParameterError(
                    f"PC-SAFT association parameter {label!r} must be finite and positive "
                    f"(got {value!r})."
                )
        for label, value in (("na", self.na), ("nb", self.nb)):
            if not np.isfinite(value) or value < 1.0:
                raise PCSAFTParameterError(
                    f"PC-SAFT association site count {label!r} must be at least 1 (got "
                    f"{value!r}); only A-B bonding is modelled, so a molecule with sites "
                    "of one type only cannot associate."
                )
        expected = self._KNOWN_SCHEMES.get(str(self.scheme).strip().upper())
        if expected is not None and (float(self.na), float(self.nb)) != expected:
            raise PCSAFTParameterError(
                f"PC-SAFT association scheme {self.scheme!r} means na, nb = {expected!r} but "
                f"the record says {(self.na, self.nb)!r}."
            )


def _resolve_segment_number(record: "PCSAFTRecord") -> float:
    """Return ``m``, derived from ``segments_per_g * MW_g_mol`` when that is how it was given.

    Raises:
        PCSAFTParameterError: If both specifications are present, if neither
            is, or if ``segments_per_g`` was given without a molar mass to
            multiply it by.
    """
    name = record.name
    if record.segments_per_g is None:
        if record.m is None:
            raise PCSAFTParameterError(
                f"PC-SAFT record for {name!r} gives neither 'm' nor 'segments_per_g'; one of "
                "the two is required (see ADR-0022 for the segments-per-mass convention)."
            )
        return float(record.m)

    if record.m is not None:
        raise PCSAFTParameterError(
            f"PC-SAFT record for {name!r} gives both 'm' and 'segments_per_g'. They are two "
            "spellings of the same parameter (m = segments_per_g * MW_g_mol); give exactly one."
        )
    per_gram = float(record.segments_per_g)
    if not np.isfinite(per_gram) or per_gram <= 0.0:
        raise PCSAFTParameterError(
            f"PC-SAFT parameter 'segments_per_g' for {name!r} must be finite and positive "
            f"(got {record.segments_per_g!r})."
        )
    if record.MW_g_mol is None:
        raise PCSAFTParameterError(
            f"PC-SAFT record for {name!r} gives 'segments_per_g' (mol/g) but no 'MW_g_mol'; the "
            "segment number of a polymer is m = segments_per_g * MW_g_mol, so the molar mass of "
            "this particular sample is required."
        )
    molar_mass = float(record.MW_g_mol)
    if not np.isfinite(molar_mass) or molar_mass <= 0.0:
        raise PCSAFTParameterError(
            f"PC-SAFT parameter 'MW_g_mol' for {name!r} must be finite and positive "
            f"(got {record.MW_g_mol!r})."
        )
    return per_gram * molar_mass


@dataclass(frozen=True)
class PCSAFTRecord:
    """Pure-component PC-SAFT parameters for one compound.

    ``name`` is the canonical (normalized) component name. ``association`` is
    ``None`` for a non-associating compound and a
    :class:`PCSAFTAssociationRecord` otherwise (ADR-0018).

    Two ways to say how long the chain is (ADR-0022)
    ------------------------------------------------
    A small molecule has a segment number: ``m`` is given directly, which is
    how every packaged record and every record written before ADR-0022 is
    spelled, and ``MW_g_mol`` stays the optional, unused annotation it always
    was (the component databank is the source of truth for molar mass).

    A **polymer** does not. Its parameters are fitted per unit *mass* and the
    chain length follows from the molar mass of the particular sample, which is
    the ``m/M`` convention of Gross & Sadowski's polymer work: give
    ``segments_per_g`` in mol/g together with ``MW_g_mol`` in g/mol, and

        m = segments_per_g * MW_g_mol

    is derived here, once, so no caller multiplies it by hand. ``MW_g_mol`` is
    then **required** and load-bearing rather than annotation - it is the one
    number that distinguishes a 16 400 g/mol polyethylene (``m = 431.32``) from
    a 53 000 g/mol one (``m = 1393.9``).

    Exactly one of ``m`` and ``segments_per_g`` may be given. Both, or neither,
    is an error rather than a silently preferred default: the two say the same
    thing in different units, and a record carrying disagreeing versions of one
    parameter is a data defect, not a choice.

    This models a **monodisperse** polymer - one chain length, one component.
    A real sample has a distribution (the sources cited in
    ``tests/fixtures/pcsaft/martini2009_polymers.json`` quote polydispersities
    of 1.14 to 2.94), and representing it as several pseudo-components is a
    different capability that this record does not provide.
    """

    name: str
    m: float | None = None
    sigma_A: float | None = None
    epsilon_k_K: float | None = None
    MW_g_mol: float | None = None
    source: str = ""
    association: PCSAFTAssociationRecord | None = None
    segments_per_g: float | None = None

    def __post_init__(self) -> None:
        if not self.name:
            raise PCSAFTParameterError("PC-SAFT component names must be non-empty.")
        # ``m`` may be given directly or derived from the segments-per-mass
        # convention; ``sigma_A`` and ``epsilon_k_K`` are only defaulted so
        # that ``m`` can precede them with a default of its own, and a record
        # without them is rejected right here.
        object.__setattr__(self, "m", _resolve_segment_number(self))
        for label, value in (
            ("m", self.m),
            ("sigma_A", self.sigma_A),
            ("epsilon_k_K", self.epsilon_k_K),
        ):
            if value is None:
                raise PCSAFTParameterError(
                    f"PC-SAFT record for {self.name!r} is missing the required parameter {label!r}."
                )
            if not np.isfinite(value) or value <= 0.0:
                raise PCSAFTParameterError(
                    f"PC-SAFT parameter {label!r} for {self.name!r} must be finite and "
                    f"positive (got {value!r})."
                )
        if self.MW_g_mol is not None and (not np.isfinite(self.MW_g_mol) or self.MW_g_mol <= 0.0):
            raise PCSAFTParameterError(
                f"PC-SAFT parameter 'MW_g_mol' for {self.name!r} must be finite and "
                f"positive when given (got {self.MW_g_mol!r})."
            )


@dataclass(frozen=True)
class PCSAFTParameters:
    """Pure-component PC-SAFT parameters keyed by normalized component name."""

    records: Mapping[str, PCSAFTRecord]

    @classmethod
    def load(cls) -> "PCSAFTParameters":
        """Return the packaged parameter set (Gross & Sadowski 2001, Table 1)."""
        return cls._from_payload(_load_payload())

    @classmethod
    def _from_payload(cls, payload: Mapping[str, object]) -> "PCSAFTParameters":
        if payload.get("schema_version") != _SCHEMA_VERSION:
            raise PCSAFTParameterError(
                "Unsupported PC-SAFT parameter schema "
                f"{payload.get('schema_version')!r}; expected {_SCHEMA_VERSION}."
            )
        if str(payload.get("model", "")).strip().casefold() != "pc-saft":
            raise PCSAFTParameterError("PC-SAFT parameter payload has an unexpected model label.")

        components = payload.get("components")
        if not isinstance(components, list):
            raise PCSAFTParameterError("PC-SAFT parameter payload must define a 'components' list.")
        return cls.from_records(components)

    @classmethod
    def from_records(
        cls, records: Iterable[PCSAFTRecord | Mapping[str, object]]
    ) -> "PCSAFTParameters":
        """Build a parameter set from user-supplied records.

        Each entry is either a :class:`PCSAFTRecord` or a mapping with the keys
        ``name``, ``sigma_A``, ``epsilon_k_K``, the optional ``MW_g_mol`` /
        ``source`` / ``association``, and **exactly one** of

        - ``m`` - the segment number, as every packaged record gives it; or
        - ``segments_per_g`` (mol/g) together with ``MW_g_mol`` (g/mol), the
          polymer convention of ADR-0022, from which ``m = segments_per_g *
          MW_g_mol`` is derived.

        Giving both, or neither, raises. Names are normalized with
        ``chemthermo.data.normalize_name``; duplicates are rejected.

        ``association`` is either a :class:`PCSAFTAssociationRecord` or a
        mapping with ``kappa_ab`` and ``epsilon_ab_k_K`` and the optional
        ``na`` / ``nb`` / ``scheme`` / ``source``; omit it (or pass ``None``)
        for a non-associating component.
        """
        store: dict[str, PCSAFTRecord] = {}
        for entry in records:
            record = _coerce_record(entry)
            if record.name in store:
                raise PCSAFTParameterError(
                    f"Duplicate PC-SAFT parameters for component {record.name!r}."
                )
            store[record.name] = record
        if not store:
            raise PCSAFTParameterError("PC-SAFT parameter set must contain at least one record.")
        return cls(records=store)

    def names(self) -> tuple[str, ...]:
        """Return the canonical component names this set covers, sorted."""
        return tuple(sorted(self.records))

    def record(self, name: str) -> PCSAFTRecord:
        """Return the record for one component, raising if it is missing."""
        canonical = normalize_name(name)
        try:
            return self.records[canonical]
        except KeyError as exc:
            raise PCSAFTParameterError(f"Missing PC-SAFT parameters for: {name}.") from exc

    def for_components(self, names: Sequence[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(m, sigma_A, epsilon_k_K)`` arrays in the requested order.

        Raises:
            PCSAFTParameterError: If ``names`` is empty or any component has no
                record.
        """
        picked = self._picked(names)
        return (
            np.array([r.m for r in picked], dtype=float),
            np.array([r.sigma_A for r in picked], dtype=float),
            np.array([r.epsilon_k_K for r in picked], dtype=float),
        )

    def association_for_components(
        self, names: Sequence[str]
    ) -> tuple[PCSAFTAssociationRecord | None, ...]:
        """Return one association record (or ``None``) per component, in order.

        An all-``None`` result means the set is non-associating, which is what
        :class:`~chemthermo.eos.PCSAFTEOS` uses to skip the association term
        entirely (ADR-0018).
        """
        return tuple(record.association for record in self._picked(names))

    def _picked(self, names: Sequence[str]) -> list[PCSAFTRecord]:
        if not names:
            raise PCSAFTParameterError("PC-SAFT requires at least one component.")

        canonical = [normalize_name(name) for name in names]
        missing = [original for original, key in zip(names, canonical) if key not in self.records]
        if missing:
            missing_list = ", ".join(str(name) for name in missing)
            raise PCSAFTParameterError(f"Missing PC-SAFT parameters for: {missing_list}.")
        return [self.records[key] for key in canonical]


def _coerce_association(
    entry: PCSAFTAssociationRecord | Mapping[str, object] | None, owner: str
) -> PCSAFTAssociationRecord | None:
    """Return an association record from a record, a mapping or ``None``."""
    if entry is None:
        return None
    if isinstance(entry, PCSAFTAssociationRecord):
        return entry
    if not isinstance(entry, Mapping):
        raise PCSAFTParameterError(
            f"PC-SAFT 'association' for {owner!r} must be an object or null."
        )
    try:
        kappa_ab = float(entry["kappa_ab"])  # type: ignore[arg-type]
        epsilon_ab_k_K = float(entry["epsilon_ab_k_K"])  # type: ignore[arg-type]
        na = float(entry.get("na", 1.0))  # type: ignore[arg-type]
        nb = float(entry.get("nb", 1.0))  # type: ignore[arg-type]
    except KeyError as exc:
        raise PCSAFTParameterError(
            f"PC-SAFT 'association' for {owner!r} is missing key {exc}."
        ) from exc
    except (TypeError, ValueError) as exc:
        raise PCSAFTParameterError(
            f"PC-SAFT 'association' for {owner!r} contains non-numeric parameters."
        ) from exc
    return PCSAFTAssociationRecord(
        kappa_ab=kappa_ab,
        epsilon_ab_k_K=epsilon_ab_k_K,
        na=na,
        nb=nb,
        scheme=str(entry.get("scheme", "2B")),
        source=str(entry.get("source", "")),
    )


def _coerce_record(entry: PCSAFTRecord | Mapping[str, object]) -> PCSAFTRecord:
    if isinstance(entry, PCSAFTRecord):
        # ``entry.m`` is already resolved (``__post_init__`` derived it from
        # ``segments_per_g`` if that is how it was written), so the copy takes
        # the segment number directly and must not re-derive it.
        return PCSAFTRecord(
            name=normalize_name(entry.name),
            m=entry.m,
            sigma_A=entry.sigma_A,
            epsilon_k_K=entry.epsilon_k_K,
            MW_g_mol=entry.MW_g_mol,
            source=entry.source,
            association=entry.association,
        )
    if not isinstance(entry, Mapping):
        raise PCSAFTParameterError("PC-SAFT component entries must be objects.")
    try:
        name = normalize_name(str(entry["name"]))
        sigma_A = float(entry["sigma_A"])  # type: ignore[arg-type]
        epsilon_k_K = float(entry["epsilon_k_K"])  # type: ignore[arg-type]
    except KeyError as exc:
        raise PCSAFTParameterError(f"PC-SAFT component entry missing key {exc}.") from exc
    except (TypeError, ValueError) as exc:
        raise PCSAFTParameterError(
            "PC-SAFT component entry contains non-numeric parameters."
        ) from exc

    raw_m = entry.get("m")
    raw_segments = entry.get("segments_per_g")
    raw_mw = entry.get("MW_g_mol")
    if raw_m is None and raw_segments is None:
        raise PCSAFTParameterError(
            f"PC-SAFT component entry for {name!r} gives neither 'm' nor 'segments_per_g'."
        )
    try:
        m = None if raw_m is None else float(raw_m)  # type: ignore[arg-type]
        segments_per_g = None if raw_segments is None else float(raw_segments)  # type: ignore[arg-type]
        mw = None if raw_mw is None else float(raw_mw)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise PCSAFTParameterError(
            "PC-SAFT component entry contains non-numeric parameters."
        ) from exc

    return PCSAFTRecord(
        name=name,
        m=m,
        sigma_A=sigma_A,
        epsilon_k_K=epsilon_k_K,
        MW_g_mol=mw,
        source=str(entry.get("source", "")),
        association=_coerce_association(
            entry.get("association"),  # type: ignore[arg-type]
            name,
        ),
        segments_per_g=segments_per_g,
    )


@lru_cache(maxsize=1)
def _load_payload() -> dict[str, object]:
    """Load the raw PC-SAFT parameter payload from package resources."""
    data_path = resources.files("chemthermo.parameters")
    for segment in _PCSAFT_RESOURCE:
        data_path = data_path / segment
    with data_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@lru_cache(maxsize=1)
def get_pcsaft_parameters() -> PCSAFTParameters:
    """Return the default (packaged) PC-SAFT parameter set."""
    return PCSAFTParameters.load()
