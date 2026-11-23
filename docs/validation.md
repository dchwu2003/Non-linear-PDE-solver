# Measured test results

These results come from a run on 2026-10-23 with **MATLAB 24.2.0.3361896 (R2024b), Update 10**, on an Ubuntu GitHub-hosted runner. The [CI run](https://github.com/dchwu2003/Non-linear-PDE-solver/actions/runs/40029637502) passed the tests and the full example, including figure generation.

Code version used for the solver and the experiments: [`3a89a08501076936d17c97f0cd1213e4c75a5bb1`](https://github.com/dchwu2003/Non-linear-PDE-solver/commit/3a89a08501076936d17c97f0cd1213e4c75a5bb1). The CSV files and figures in this repository are copies of that run's output, not hand-made examples. Running the commands in the README writes new outputs to `results/` and does not overwrite these copies.

## Grid refinement at fixed ε

Both test problems use ε = 0.2 on `[0,2] x [0,1]`. Both grid spacings halve from one grid to the next. The error is the largest difference from the exact solution at the nodes. The observed order is `log2(error_coarse/error_fine)`.

| Grid | Error, x-only test | Order | Error, x–y test | Order |
|---|---:|---:|---:|---:|
| 21 x 13 | 3.776507e-4 | — | 1.656320e-4 | — |
| 41 x 25 | 9.488772e-5 | 1.9928 | 4.130739e-5 | 2.0035 |
| 81 x 49 | 2.369195e-5 | 2.0018 | 1.032057e-5 | 2.0009 |

Each x-only solve took 48 iterations, and each x–y solve took 35. Residuals recomputed directly from the finite-difference formulas were below `7.1e-12` in every case. For the x–y test, the separately estimated boundary-derivative error went from `4.2028e-4` to `6.6035e-6` over these grids. This fast decrease is specific to this smooth test function and should not be read as a general accuracy result.

![Error versus grid size](figures/grid_convergence.png)

Full-precision data: [analytic.csv](data/analytic.csv) (x-only test), [manufactured.csv](data/manufactured.csv) (x–y test). [numerical-method.md](numerical-method.md) gives the formulas and boundary conventions.

## Default example: lowering ε step by step

The example uses 121 x 61 nodes on `[0,6] x [0,5]`, left/right values 1 and -1.1, zero y derivatives on the top and bottom, relaxation 0.5, relative update tolerance `1e-8`, and absolute residual tolerance `1e-7`. Every ε value is solved. The change in the solution from one ε to the next is recorded but is not used to stop.

These settings differ from the run shown on the December 2025 poster (boundary values +1/−1, ε down to 0.01, 320×270 grid).

| ε | Iterations | Final residual | Maximum cell Péclet number |
|---|---:|---:|---:|
| 0.1000 | 649 | 9.8386e-8 | 0.275 |
| 0.0875 | 263 | 9.5291e-8 | 0.314 |
| 0.0750 | 269 | 9.8715e-8 | 0.367 |
| 0.0625 | 277 | 9.5936e-8 | 0.440 |
| 0.0500 | 287 | 9.5548e-8 | 0.550 |

![Residual history for each ε](figures/residuals.png)

At the first ε, the residual stays nearly flat for a long time before it drops. A low iteration limit would stop this example too early. The sparse matrices are built once per ε, but the iteration is still fairly slow here. We make no claims about speed.

![Solution profiles at mid-height](figures/profiles.png)

![Final solution surface](figures/solution.png)

Because the right boundary value is more negative than the left value is positive, the transition sits near the left boundary. For these boundary data the solution does not depend on y. The small residuals show that the iteration solved the discrete equations. They do not measure how accurate the grid is in the thin layer. The grid refinement results above apply only to those two test problems and their settings.

Full-precision data: [continuation.csv](data/continuation.csv). The CI artifact `matlab-validation` contains the full solution arrays and histories in `example.mat` and the test report in `validation.mat`, subject to GitHub's artifact retention policy. The figures and CSV copies stay in the repository.

## Other checks

The same run also checked:

- zero and nonzero constant solutions;
- a linear field with nonzero top/bottom derivatives and a perturbed starting guess on a non-square grid;
- residuals recomputed from the formulas, and errors in the fixed boundary values;
- that the solver does not report convergence with a tiny relaxation and a large residual;
- rejection of invalid ε, boundary lengths and relaxation inputs;
- that the boundary-derivative error of the x–y test decreases under refinement.

CI reports success only if all assertions pass. The workflow runs on pull requests and branch pushes, using MATLAB R2024b without extra toolboxes.
