function [U_final, info] = Nonlinear_iteration(U_initial, epsilon, h, p, Nx, Ny, ...
    bottom, top, left, right, max_iterations, error_tolerance, options)
%NONLINEAR_ITERATION Damped Picard solve of U*Ux-epsilon*Uxx-Uyy=f.
% Existing 12-argument calls remain supported. Optional fields:
% forcing (Nx-by-Ny), relaxation (0,1], residual_tolerance (absolute infinity
% norm on PDE rows), verbose. All vectors use x-fast ordering throughout.
% Nonconvergence returns info.converged=false and emits a warning.
    if nargin < 13, options = struct(); end
    defaults = struct('forcing', zeros(Nx,Ny), 'relaxation', 0.5, ...
        'residual_tolerance', error_tolerance, 'verbose', true);
    names = fieldnames(options);
    for k = 1:numel(names)
        if ~isfield(defaults, names{k}), error('pde:UnknownOption','Unknown option: %s',names{k}); end
        defaults.(names{k}) = options.(names{k});
    end
    options = defaults;
    validateattributes(max_iterations, {'numeric'}, {'scalar','integer','positive','finite'});
    validateattributes(error_tolerance, {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(options.relaxation, {'numeric'}, {'scalar','real','finite','>',0,'<=',1});
    validateattributes(options.residual_tolerance, {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(options.verbose, {'logical','numeric'}, {'scalar','binary'});
    validateattributes(U_initial, {'numeric'}, {'vector','numel',Nx*Ny,'real','finite'});
    [D, b, Dx, free] = assemble_system(epsilon,h,p,Nx,Ny,bottom,top,left,right,options.forcing);
    u = U_initial(:);
    ids = reshape(1:Nx*Ny,Nx,Ny);
    u(ids(1,:)) = left(:); u(ids(end,:)) = right(:);
    info.converged = false;
    info.relative_update = zeros(max_iterations,1);
    info.residual = zeros(max_iterations,1);
    info.initial_residual = norm((D*u + u.*(Dx*u) - b),inf);
    for iteration = 1:max_iterations
        % Dx has zero Dirichlet rows, so their identity equations are preserved.
        A = D + spdiags(u,0,Nx*Ny,Nx*Ny)*Dx;
        candidate = A\b;
        next = u + options.relaxation*(candidate-u);
        if any(~isfinite(next)), error('pde:NonfiniteIterate','Nonfinite linear solve.'); end
        update = norm(next(free)-u(free),inf)/max(1,norm(next(free),inf));
        residual = D*next + next.*(Dx*next) - b;
        info.relative_update(iteration) = update;
        info.residual(iteration) = norm(residual(free),inf);
        u = next;
        if options.verbose
            fprintf('Picard %d: relative update %.3e, PDE residual %.3e\n', ...
                iteration, update, info.residual(iteration));
        end
        if update <= error_tolerance && info.residual(iteration) <= options.residual_tolerance
            info.converged = true;
            break;
        end
    end
    U_final = u;
    info.iterations = iteration;
    info.relative_update = info.relative_update(1:iteration);
    info.residual = info.residual(1:iteration);
    info.dirichlet_error = max(abs([u(ids(1,:))-left(:); u(ids(end,:))-right(:)]));
    M = reshape(u,Nx,Ny);
    % Independent one-sided derivative estimates; these include truncation error
    % and need not be zero even though the eliminated ghost equations are exact.
    gb = (-3*M(:,1)+4*M(:,2)-M(:,3))/(2*p);
    gt = (3*M(:,end)-4*M(:,end-1)+M(:,end-2))/(2*p);
    info.neumann_error = max(abs([gb(2:end-1)-reshape(bottom(2:end-1),[],1); ...
        gt(2:end-1)-reshape(top(2:end-1),[],1)]));
    info.max_cell_peclet = max(abs(u(free)))*h/(2*epsilon);
    if ~info.converged
        warning('pde:NotConverged', 'Picard stopped after %d iterations; PDE residual %.3e.', ...
            iteration, info.residual(end));
    end
end
