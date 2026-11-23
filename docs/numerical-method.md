# Numerical method and tests

This page gives the formulas behind the README. Where a standard technical name exists, it is given after a plain description, for readers who want to look it up.

We solve

$$u u_x-\varepsilon u_{xx}-u_{yy}=f.$$

On the left and right edges, the value of $u$ is fixed (Dirichlet conditions). On the bottom and top edges, the derivative of $u$ in the **positive y direction** is fixed (Neumann conditions). The default example uses $f=0$. The optional source term $f$ is only there for test problems: you choose an exact solution, compute the $f$ it needs, and check that the solver recovers it (the method of manufactured solutions).

## Grid layout and operator

`U` is an `Nx`-by-`Ny` array, with x varying along its first dimension. As a vector it is always `u=U(:)`, with global index `i+(j-1)*Nx`. For plotting, use `U'`. The iteration never changes this layout.

Interior derivatives use centered differences. For the bottom derivative $g_b$, we remove a ghost value outside the grid using

$$u_{i,0}=u_{i,2}-2h_y g_b(x_i).$$

Here the physical bottom boundary has index 1. The resulting bottom row has coefficient $-2/h_y^2$ for its upper neighbor and right-hand-side contribution $-2g_b/h_y$. At the top, the corresponding right-hand-side contribution is $+2g_t/h_y$. The four corner equations use the fixed left/right values, so the derivative conditions apply to boundary points other than the corners.

`assemble_system` builds the second-derivative matrix `D`, the x-derivative matrix `Dx`, and the right-hand side once for each ε. Rows of `Dx` for fixed-value boundary points are zero. Each iteration freezes the current `u` in the nonlinear term and solves the linear system

```
(D + diag(u_old)*Dx)*candidate = b
u_new = u_old + relaxation*(candidate-u_old)
```

with a sparse direct solve. The update moves only part of the way toward `candidate`. The default fraction is 0.5. In standard terms, this is a damped (relaxed) Picard iteration.

Lowering ε step by step and starting each step from the previous solution is called parameter continuation. The grid stays the same throughout. This is not a multigrid method.

## Stopping rule and diagnostics

The iteration stops only when both conditions hold:

- `norm(u_new-u_old,inf)/max(1,norm(u_new,inf)) <= update_tolerance` on points that are not fixed-value boundary points;
- `norm(D*u_new + u_new.*(Dx*u_new)-b,inf) <= residual_tolerance` on equation rows.

The residual tolerance is absolute, so its meaning depends on the scaling of the equation. `info.dirichlet_error` checks the fixed boundary values separately. `info.neumann_error` estimates the boundary derivatives separately with second-order one-sided differences. It includes truncation error, so it should decrease as the grid is refined, but it will not drop to the solver tolerance. If the maximum number of iterations is reached, the code raises `pde:NotConverged`, and the example treats this as an error.

`info.max_cell_peclet = max(abs(u))*h_x/(2*epsilon)` on equation rows compares the carrying term with the smoothing term over one grid cell (the cell Péclet number). Values above one can make the centered-difference coefficients change sign and produce wiggles. A value below one is a useful check on the grid, but it does not prove accuracy or convergence. We do not claim the method always converges or stays accurate for small ε.

## Test problems

`tests/run_tests.m` uses plain assertions. It does not need MATLAB's unit-testing toolbox.

1. Constant solutions, including all-zero boundary values.
2. A linear field `u=0.2+0.1*x+0.05*y`, with nonzero top/bottom derivatives and `f=0.1*u`, on a non-square grid. The starting guess is perturbed, so the test runs several nonlinear updates and checks the vector ordering.
3. A very small relaxation with a single iteration. A small update must not count as converged when the residual is still large.
4. Invalid ε, boundary lengths, and relaxation inputs.
5. A known nonlinear solution, at fixed ε = 0.2:

   $$u(x,y)=-0.5\tanh\left(\frac{0.5(x-1)}{2\varepsilon}\right),\quad (x,y)\in[0,2]\times[0,1].$$

   Its values at the ends give the left/right boundary values. Its y derivatives are zero.

6. A test solution that depends on both x and y, on the same domain:

   $$u=0.2+0.1x+0.05y+0.02y^2+0.03\sin(\pi x/2)\cos(\pi y).$$

   Its exact derivatives define `f=u*ux-epsilon*uxx-uyy`, with bottom derivative 0.05 and top derivative 0.09. This tests both signs of the boundary-derivative term and a solution that varies in y.

Both refinement studies use grids of 21-by-13, 41-by-25, and 81-by-49 nodes, halving both spacings each time. The error is the largest absolute difference from the exact solution at the nodes. The observed order is `log2(error_coarse/error_fine)`. The tests require orders between 1.8 and 2.2. They also recompute each solution's residual from the finite-difference formulas directly, without using the assembled matrix.

These tests only show how the code behaves on smooth test solutions. They do not show second-order accuracy in general, for example for unresolved thin layers or other boundary data. Changing ε changes the equation itself, so it is kept separate from the grid refinement checks.

## Reproducing the outputs

From the repository root:

```matlab
addpath('tests');
report = run_tests('results', true);
result = run_example('results', true);
```

To run without graphics, set the second argument to `false`. The CSV files contain errors, observed orders, residuals and iteration counts. The MAT files keep the settings and arrays. GitHub Actions uploads these files and figures as the `matlab-validation` artifact. [validation.md](validation.md) records the run details and results for the reviewed version.
