# Nonlinear PDE Solver in MATLAB

This project solves the nonlinear partial differential equation

\[
U \, U_x - U_{yy} - \varepsilon U_{xx} = 0, \quad (x,y) \in (0,20)\times(0,15)
\]

subject to mixed boundary conditions:
- \( U(0,y) = 1 \)
- \( U(20,y) = -1 \)
- \( U_y(x,0) = U_y(x,15) = 0 \)

The solver employs a **finite difference discretization** and **nonlinear iteration scheme**, with multilevel refinement and variable diffusion parameter \(\varepsilon\).

---

## Features

- Multi-resolution mesh refinement (4 levels by default)  
- Nonlinear iteration for \( U U_x \) term  
- Mixed Dirichlet–Neumann boundary handling  
- Automatic prolongation of coarse-grid solutions  
- Visualization via surface plots

---

## File Overview

| File | Description |
|------|--------------|
| `main.m` | Main script that runs the PDE solver on multiple grid levels |
| `Initial_guess.m` | Builds the linearized system for the first iteration |
| `nonlinear_iteration.m` | Performs iterative updates until convergence |

---

## How to Run

Open MATLAB and execute:
```matlab
run('main.m')
```
You should see a series of 3D surface plots, one per refinement level.

Expected console output:
```
Precomputing initial guess with eps = 0.20 on (40 x 30) grid...
===== Level 1 =====
...
Nonlinear iteration stops at step 9, rel err = 9.1e-05
```

---

## Parameters

| Parameter | Meaning | Default |
|------------|----------|----------|
| `Nx0`, `Ny0` | Initial grid size | 100, 100 |
| `eps0` | Initial diffusion coefficient | 0.1 |
| `n_levels` | Number of refinement levels | 4 |
| `max_iter` | Maximum nonlinear iterations | 20 |
| `err_ctrl` | Convergence tolerance | 1e-4 |

---

## Dependencies

- MATLAB R2021a or later  
- No external toolboxes required

---

## Citation / Reuse

If you adapt or extend this code, please credit the author and link to this repository.  
This solver can serve as a base for studying **nonlinear convection–diffusion PDEs** or for extending to 2D steady Navier–Stokes–type systems.
Contributors are Daocheng Wu, Shen Chang, Jiayi Tang, Yu wang, Yilin Chen

```
@software{nonlinear_pde_solver,
  author = {Daocheng Wu, Shen Chang, Jiayi Tang, Yu wang, Yilin Chen},
  title = {Nonlinear PDE Solver in MATLAB},
  year = {2025},
  url = {https://gitlab.umich.edu/logm/Fall25/nonlinear-pde-final-code/matlab-code-for-pde.git},
  note = {Finite difference nonlinear PDE solver with multigrid refinement.}
}
```

---

## License

MIT License — free to modify and redistribute with attribution.

