# Discretization and validation

We solve

$$u u_x-\varepsilon u_{xx}-u_{yy}=f,$$

with Dirichlet data at the left/right boundaries and prescribed **positive-y derivatives** at the bottom/top boundaries. The demonstration uses $f=0$. The optional forcing term exists to verify the solver with manufactured solutions.

## Layout and operator

`U` is an `Nx`-by-`Ny` array, with x varying along its first dimension. The sole vector representation is `u=U(:)`, with global index `i+(j-1)*Nx`. A plotting array is `U'`. No iteration changes this layout.

Interior derivatives use centered differences. For bottom derivative $g_b$, eliminating a ghost value gives

$$u_{i,0}=u_{i,2}-2h_y g_b(x_i).$$

Here the physical bottom boundary has index 1. The resulting bottom diffusion row has coefficient $-2/h_y^2$ for its northern neighbor and RHS contribution $-2g_b/h_y$. At the top, the corresponding RHS contribution is $+2g_t/h_y$. Dirichlet equations replace the four corner equations, so Neumann data are enforced on boundary points excluding corners.

`assemble_system` builds the diffusion matrix `D`, derivative `Dx`, and RHS once per epsilon. Dirichlet rows of `Dx` are zero. Picard iteration solves

```
(D + diag(u_old)*Dx)*candidate = b
u_new = u_old + relaxation*(candidate-u_old)
```

using a sparse direct solve. The default relaxation is 0.5. This is fixed-grid parameter continuation, not a multigrid algorithm.

## Stopping and diagnostics

Both conditions must hold:

- `norm(u_new-u_old,inf)/max(1,norm(u_new,inf)) <= update_tolerance` on non-Dirichlet nodes;
- `norm(D*u_new + u_new.*(Dx*u_new)-b,inf) <= residual_tolerance` on PDE rows.

The residual tolerance is absolute and has the units/scaling of the chosen equation. `info.dirichlet_error` checks boundary values separately. `info.neumann_error` uses independent second-order one-sided derivative estimates; it includes truncation error and is expected to decrease with refinement, rather than vanish at solver tolerance. Maximum iterations produce `pde:NotConverged`, and the example treats a failed stage as an error.

`info.max_cell_peclet = max(abs(u))*h_x/(2*epsilon)` on PDE rows indicates convection resolution. Values above one can destroy the sign structure of centered-difference coefficients. A value below one is a useful stencil check, not proof of accuracy or convergence. No unconditional convergence or small-epsilon accuracy is claimed.

## Verification problems

`tests/run_tests.m` uses plain assertions, independent of MATLAB's unit-testing toolbox.

1. Constant solutions, including zero boundary values.
2. An affine field `u=0.2+0.1*x+0.05*y`, with nonzero Neumann data and `f=0.1*u`, on a rectangular grid. A perturbed initial iterate exercises multiple nonlinear updates and vector ordering.
3. A deliberately tiny relaxation and one iteration: a small update must not pass when the PDE residual remains large.
4. Invalid epsilon, boundary length, and relaxation inputs.
5. A nonlinear analytic solution, at fixed epsilon 0.2:

   $$u(x,y)=-0.5\tanh\left(\frac{0.5(x-1)}{2\varepsilon}\right),\quad (x,y)\in[0,2]\times[0,1].$$

   Its endpoint values provide Dirichlet data and its y derivatives vanish.

6. A genuinely two-dimensional manufactured solution on the same domain:

   $$u=0.2+0.1x+0.05y+0.02y^2+0.03\sin(\pi x/2)\cos(\pi y).$$

   Analytic derivatives define `f=u*ux-epsilon*uxx-uyy`, with bottom derivative 0.05 and top derivative 0.09. This exercises both signs of the Neumann contribution and y-dependent solution values.

Both refinement studies use 21-by-13, 41-by-25, and 81-by-49 nodes, halving both spacings. Reported error is the maximum absolute nodal error. Observed order is `log2(error_coarse/error_fine)`. The suite requires orders between 1.8 and 2.2 and independently evaluates explicit PDE stencils, without using the assembled matrix, to verify each solution's discrete residual.

These tests demonstrate behavior for smooth benchmark solutions. They do not establish universal second-order convergence for unresolved internal layers or all boundary data. Epsilon continuation changes the PDE and is deliberately kept separate from grid convergence.

## Reproducing outputs

From the repository root:

```matlab
addpath('tests');
report = run_tests('results', true);
result = run_example('results', true);
```

The same calculations can run without graphics by setting the second argument to `false`. CSVs contain errors, observed orders, residuals and iteration counts; MAT files retain parameters and arrays. GitHub Actions uploads these files and figures as the `matlab-validation` artifact. Runtime and results for the reviewed version are recorded in [validation.md](validation.md).
