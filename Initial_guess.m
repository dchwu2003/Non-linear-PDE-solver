function U_initial = Initial_guess(epsilon, h, p, Nx, Ny, bottom, top, left, right)
%INITIAL_GUESS Solve -epsilon*Uxx-Uyy=0 with the solver's boundary stencil.
% Returns an x-fast vector: reshape(U_initial,Nx,Ny) has x along rows.
    [A, b] = assemble_system(epsilon, h, p, Nx, Ny, bottom, top, left, right);
    U_initial = A\b;
end
