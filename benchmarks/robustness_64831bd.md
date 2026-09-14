# Robustness map at `64831bd`

- Generated: 2026-09-14T16:12:50+00:00
- Scope: full sweep, all families
- Machine: Apple M2 Max / Python 3.11.6 / numpy 2.4.2
- Total wall time: 2202.5 s
- Tree dirty at measurement: True

What this measures, and what it does not: ADR-0027.

## Family x outcome

| family | states | verdicts | refusal classes | invariant violations | worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `pr-phi-phi` | 1270 | VLE 319, single-liquid 406, single-vapor 545 | - | 0 | 2.12e-13 | 7.79e-08 | -1.20e-05 | 5.5 |
| `pcsaft` | 204 | VLE 132, single-liquid 72 | - | 0 | 1.86e-13 | 3.58e-08 | -3.95e-04 | 70.5 |
| `pcsaft-associating` | 260 | LLE 145, VLE 5, single-liquid 89, single-vapor 21 | - | 0 | 5.07e-13 | 2.83e-09 | -5.34e-07 | 456.5 |
| `modified-raoult` | 108 | LLE 48, VLE 28, VLLE 4, single-liquid 13, single-vapor 15 | - | 0 | 1.11e-16 | 9.23e-13 | -1.07e-06 | 3.7 |
| `gamma-gamma` | 16 | LLE 15, single-liquid 1 | - | 0 | 2.22e-16 | 1.43e-13 | -1.07e-06 | 0.7 |
| `polymer` | 252 | LLE 140, VLE 50, single-liquid 62 | - | 0 | 2.22e-16 | 2.09e-07 | -4.12e-07 | 346.3 |
| `eos-three-phase` | 117 | LLE 46, LLL 12, VLE 27, VLLE 19, single-liquid 10, single-vapor 3 | - | 0 | 2.11e-12 | 2.80e-08 | -1.25e-04 | 655.8 |
| `gamma-phi-legacy` | 30 | VLE 8, single-liquid 3, single-vapor 14 | rr-no-bracket 5 | 0 | 1.02e-13 | 0.00e+00 | - | 0.0 |
| `pr-near-critical` | 104 | VLE 55, single-liquid 28, single-vapor 21 | - | 0 | 1.47e-13 | 3.33e-08 | -1.94e-07 | 2.0 |
| `pcsaft-associating-ternary` | 144 | LLE 84, VLE 8, VLLE 6, single-liquid 30, single-vapor 16 | - | 0 | 1.74e-13 | 2.98e-09 | -1.30e-04 | 661.5 |

**Totals:** 2505 states, 2500 converged, 5 refused, 0 converged-but-violating, 0 skipped.

## Refusal classes, ranked by count

| refusal class | count | families | example state | message (first line) |
| --- | ---: | --- | --- | --- |
| `rr-no-bracket` | 5 | gamma-phi-legacy | `Methane/Ethane` z=(0.5, 0.5) T=200 K P=3e+06 Pa | Rachford-Rice failed to bracket a vapor fraction. |

### ... and by stage, which is what names the defect

| refusal class | stage | count | system | example state |
| --- | --- | ---: | --- | --- |
| `rr-no-bracket` | `phi-phi` | 5 | `gammaphi-methane-ethane` | `Methane/Ethane` z=(0.5, 0.5) T=200 K P=3e+06 Pa |

## Slowest states

| system | T / K | P / Pa | outcome | time / s |
| --- | ---: | ---: | --- | ---: |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 333 | 1.01e+05 | VLLE | 40.40 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 333 | 1.01e+05 | VLLE | 35.99 |
| `assoc3-water-ethanol-n-hexane` | 340 | 1e+05 | VLLE | 35.90 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 333 | 1.01e+05 | VLLE | 29.91 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 335 | 1.01e+05 | VLLE | 28.32 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 337 | 1.01e+05 | VLLE | 28.16 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 335 | 1.01e+05 | VLLE | 25.90 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 333 | 1.01e+05 | VLLE | 25.35 |
| `assoc3-water-ethanol-n-hexane` | 340 | 1e+05 | VLLE | 24.87 |
| `eos3p-pcsaft-water-ethanol-n-hexane` | 337 | 1.01e+05 | VLLE | 24.23 |
