# Four-bed SO2 converter

Catalyst bed sizing and inventory optimisation for the converter of a
150 MTPD double-contact sulfuric acid plant. Four adiabatic beds over
V2O5 on silica gel, with interbed cooling and interpass absorption of SO3
between beds 3 and 4.

## Model

Each bed is an adiabatic, pseudo-homogeneous, isobaric plug flow reactor.
The SO2 and energy balances over a slice `dz` give

```
(CA0 u0 / rhoB) df/dz = rm
G cp dT               = rm rhoB (-dH) dz
```

Eliminating `dz` and integrating with `dH(T) = -dHa + dHb T`:

```
T = (T0 + a (f - f0)) / (1 + b (f - f0))
```

valid in kelvin only, since `dH` is linear in absolute temperature. Depth
follows from `dz = alpha df`, `alpha = CA0 u0 / (rhoB rm)`. The rate keeps
the reverse step, so `alpha` diverges at a finite ceiling conversion beyond
which no bed depth is sufficient.

## Files

| file | contents |
|---|---|
| `ConverterModel.m` | basis, kinetics, bed march, ceiling, pressure drop |
| `ConverterOptimization.m` | minimum catalyst at fixed conversion |
| `run_converter.m` | driver: sizing, then optimisation |
| `test_converter.m` | regression and invariant tests |

## Usage

```matlab
run_converter                          % sizing and optimisation
runtests('test_converter.m')           % tests
```

MATLAB R2016b or later. Sizing needs no toolboxes; the optimisation needs
`fmincon` (Optimization Toolbox).

## Result

| | design | optimum |
|---|---|---|
| bed depths [m] | 0.513 / 1.214 / 1.298 / 1.499 | 0.500 / 0.500 / 1.022 / 1.648 |
| total [m] | 4.524 | 3.670 |
| catalyst [kg] | 17 588 | 14 268 |
| conversion | 0.99875 | 0.99875 |
| peak temperature [C] | 616.4 | 619.6 |

## Assumptions worth knowing

Catalyst effectiveness factor is 1 and `cp` is constant, so the beds are
sized on an apparent rate at a fixed particle size. Beds 1 to 3 run above
99 per cent of their ceiling conversion, where depth is sensitive to the
kinetic parameters and should not be read to four figures. The optimum is
held by the minimum bed depth `b.Lmin`, which is a judgement about flow
distribution rather than a derived quantity.
