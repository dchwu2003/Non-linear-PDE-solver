function config = default_config()
%DEFAULT_CONFIG Reproducible demonstration for the nonlinear equation.
    config.Nx = 320;
    config.Ny = 270;
    config.Lx = 6;
    config.Ly = 5;
    config.epsilon = linspace(0.1,0.01,10);
    config.left = 1;
    config.right = -1;
    config.max_iterations = 1000;
    config.update_tolerance = 1e-8;
    config.residual_tolerance = 1e-7;
    config.relaxation = 0.5;
end
