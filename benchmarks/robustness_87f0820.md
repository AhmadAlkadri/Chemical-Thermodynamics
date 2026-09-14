# Robustness map at `87f0820`

- Generated: 2026-09-14T09:54:49+00:00
- Scope: full sweep, all families
- Machine: Apple M2 Max / Python 3.11.6 / numpy 2.4.2
- Total wall time: 865.6 s
- Tree dirty at measurement: True

What this measures, and what it does not: ADR-0027.

## Family x outcome

| family | states | verdicts | refusal classes | invariant violations | worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `pr-phi-phi` | 1270 | VLE 319, single-liquid 406, single-vapor 545 | - | 0 | 2.12e-13 | 7.79e-08 | -1.20e-05 | 5.6 |
| `pcsaft` | 204 | VLE 132, single-liquid 72 | - | 0 | 1.86e-13 | 3.58e-08 | -3.95e-04 | 72.7 |
| `pcsaft-associating` | 260 | LLE 145, VLE 5, single-liquid 89, single-vapor 21 | - | 0 | 5.07e-13 | 2.83e-09 | -5.34e-07 | 470.9 |
| `modified-raoult` | 108 | LLE 48, VLE 28, VLLE 4, single-liquid 13, single-vapor 15 | - | 0 | 1.11e-16 | 9.23e-13 | -1.07e-06 | 3.8 |
| `gamma-gamma` | 16 | LLE 15, single-liquid 1 | - | 0 | 2.22e-16 | 1.43e-13 | -1.07e-06 | 0.7 |
| `polymer` | 252 | LLE 113, VLE 43, single-liquid 60 | split-non-convergence 34, stability-inconclusive 2 | 0 | 2.22e-16 | 2.09e-07 | -4.12e-07 | 311.9 |

**Totals:** 2110 states, 2074 converged, 36 refused, 0 converged-but-violating, 0 skipped.

## Refusal classes, ranked by count

| refusal class | count | families | example state | message (first line) |
| --- | ---: | --- | --- | --- |
| `split-non-convergence` | 34 | polymer | `Polyethylene/n-Pentane` z=(1.375e-05, 1) T=453 K P=3.6e+06 Pa | flash_tp did not converge the phi-phi split; equal-fugacity residual=3.897e-05 after 100 successive-substitution (max_delta_k=1.539e+89) and 101 second-order iterations. |
| `stability-inconclusive` | 2 | polymer | `Polyethylene/n-Pentane` z=(0.0002402, 0.9998) T=453 K P=1.05e+07 Pa | Tangent-plane stability analysis was inconclusive (no trial converged), so flash_tp cannot decide whether the feed is one phase or two |

### ... and by stage, which is what names the defect

| refusal class | stage | count | system | example state |
| --- | --- | ---: | --- | --- |
| `split-non-convergence` | `phi-phi` | 26 | `polymer-pe53000-pentane` | `Polyethylene/n-Pentane` z=(1.375e-05, 1) T=453 K P=3.6e+06 Pa |
| `split-non-convergence` | `beta-outside-window` | 4 | `polymer-pe16400-pentane` | `Polyethylene/n-Pentane` z=(4.443e-05, 1) T=453 K P=3e+05 Pa |
| `split-non-convergence` | `log-space` | 3 | `polymer-pe53000-pentane` | `Polyethylene/n-Pentane` z=(1.375e-05, 1) T=453 K P=3e+05 Pa |
| `stability-inconclusive` | `feed` | 2 | `polymer-pe53000-pentane` | `Polyethylene/n-Pentane` z=(0.0002402, 0.9998) T=453 K P=1.05e+07 Pa |
| `split-non-convergence` | `phi-phi` | 1 | `polymer-pe16400-pentane` | `Polyethylene/n-Pentane` z=(0.0002315, 0.9998) T=453 K P=7.5e+06 Pa |

## Slowest states

| system | T / K | P / Pa | outcome | time / s |
| --- | ---: | ---: | --- | ---: |
| `assoc-water-1-propanol-n-hexane` | 290 | 1e+05 | LLE | 17.17 |
| `polymer-pe53000-pentane-hexane` | 453 | 1e+06 | VLE | 15.05 |
| `assoc-water-ethanol` | 350 | 1e+05 | LLE | 12.94 |
| `assoc-water-1-propanol-n-hexane` | 290 | 5e+05 | LLE | 10.91 |
| `assoc-water-ethanol` | 350 | 1e+06 | LLE | 7.73 |
| `assoc-water-ethanol` | 350 | 5e+05 | LLE | 7.50 |
| `assoc-water-ethanol` | 350 | 2e+06 | LLE | 7.10 |
| `polymer-pe16400-pentane-hexane` | 453 | 1e+06 | VLE | 6.70 |
| `polymer-pe53000-pentane-hexane` | 453 | 5e+06 | LLE | 6.54 |
| `polymer-pe53000-pentane-hexane` | 453 | 3e+06 | LLE | 6.46 |
