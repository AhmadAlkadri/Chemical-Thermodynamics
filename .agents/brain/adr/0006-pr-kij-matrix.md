# ADR-0006: Per-pair `kij` for `PengRobinsonEOS`, and a diagonal bug fix

Status: accepted
Date: 2026-09-13

## Context
`PengRobinsonEOS` was a frozen dataclass with `kij: float = 0.0`. Both
`fugacity_coefficients` and `compressibility_factor` computed

```python
aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - self.kij)
```

which scales the DIAGONAL of `aij` (the pure-component energy parameter
`a_ii`) by `(1 - kij)` as well as every off-diagonal (cross) term. With any
nonzero `kij` this silently corrupts the pure-component limit: a binary
mixture's `a_mix` at `y -> [1, 0]` no longer equals `a_11`, so mixture and
pure-component fugacities become thermodynamically inconsistent, and any
downstream consumer (`flash_tp`, `stability_tp`) inherits the error. There was
also no way to give different `kij` values to different component pairs; a
single scalar applied uniformly to the whole mixture. ADR-0005 recorded this
as a known limitation and deferred it to this slice.

## Decision
1. `PengRobinsonEOS.kij` accepts either:
   - a `float`, applied to every off-diagonal pair `i != j` (never to the
     diagonal, whose `(1 - kij)` factor is always `1`); or
   - a `Mapping[tuple[str, str], float]` from an unordered pair of component
     names to a per-pair value, e.g. `{("Methane", "n-Decane"): 0.0411}`.
     Names are normalized with `chemthermo.data.normalize_name`. Both orders
     of a pair name the same value; giving both orders with *different*
     values raises `ModelError`; a pair naming the same component twice
     raises `ModelError`. A pair naming components absent from a given
     mixture is not an error -- the EOS does not know the mixture at
     construction time, so it is simply never looked up; a mixture pair with
     no matching entry defaults to `0.0`.
   - The mapping form is normalized in `__post_init__` (via
     `object.__setattr__`, since the dataclass stays frozen) into a sorted
     tuple of `((name_a, name_b), value)` entries, so equality, hashing and
     repr stay deterministic.
2. Two private helpers replace the previously duplicated mixing-rule code in
   `fugacity_coefficients` and `compressibility_factor`:
   - `_kij_matrix(mixture)` builds the dense `n x n` kij matrix in the
     mixture's own component order (zero diagonal always).
   - `_mixture_parameters(mixture, temperature, y)` returns
     `(a_i, b_i, aij, a_mix, b_mix)` from that matrix. Public method
     signatures are unchanged.
3. The rounded Peng-Robinson constants `0.45724` / `0.07780` are kept exactly
   as they are (ADR/decision already recorded; a separate slice may revisit
   this). This slice does not change them.

## Alternatives considered
- **A positional matrix aligned to mixture order**
  (`kij: Sequence[Sequence[float]]`, indexed by the order components appear
  in the `Mixture`). Rejected: fragile to reordering -- the same physical kij
  data silently means something different if a caller lists components in a
  different order, and there is no way to validate at construction time that
  the matrix "means" what the caller intended. A name-keyed mapping is
  permutation-invariant by construction, which this slice's tests exercise
  directly (`tests/test_pr_eos.py::test_pr_eos_kij_permutation_invariance`,
  `tests/validation/test_pr_kij_vs_thermo.py::test_pr_eos_kij_permutation_invariance_with_matrix_kij`).
- **A separate parameters registry object**, analogous to `NRTLParameters`
  (packaged JSON, `for_mixture`/`for_components` lookup). Deferred: that
  pattern is worth it once there is a packaged kij dataset with real
  provenance to justify a registry and a schema version; today there is no
  such packaged dataset, and a registry over zero real entries is
  scaffolding. The per-pair mapping keeps the same name-keyed, order-agnostic
  shape, so migrating to a registry later is additive, not a breaking change.

## Consequences
- Positive: pure-component behavior is now independent of `kij` by
  construction (the diagonal factor is always `1`), for any scalar or mapping
  value. Verified directly in
  `tests/test_pr_eos.py::test_pr_eos_kij_never_corrupts_the_diagonal` and
  cross-checked against `thermo`'s PRMIX in
  `tests/validation/test_pr_kij_vs_thermo.py::test_pr_eos_kij_matches_thermo_at_pure_component_limit`.
- Positive: `stability_tp` and `flash_tp` results at nonzero `kij` become
  trustworthy for the first time (both consume `PengRobinsonEOS` only through
  the existing `fugacity_coefficients`/`compressibility_factor` interface, so
  neither module needed to change).
- Breaking behavior change: a caller who was already using a nonzero scalar
  `kij` gets different (now-correct) numbers after this change. `kij = 0.0`
  (the default) is bit-identical to before
  (`tests/test_pr_eos.py::test_pr_eos_default_kij_matches_explicit_zero`,
  `tests/test_pr_eos.py::test_pr_eos_flash_two_phase_split_regression_unchanged`).
- Tradeoff: `PengRobinsonEOS.kij`'s declared type is a 3-shape union (`float`
  input, `Mapping` input, and the canonical `KijPairs` tuple the field holds
  after `__post_init__`), which is documented inline in
  `src/chemthermo/models/peng_robinson.py` rather than hidden.

## Supersedes (optional)
None. Removes the known-limitation note in ADR-0005's Consequences /
"What future slices will change" sections (superseded by this decision, not a
formal supersession).

## Superseded by (optional)
None.
