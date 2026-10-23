function report = run_tests(output_dir, make_plots)
%RUN_TESTS Numerical regression suite, runnable in MATLAB or GNU Octave.
% From repository root: addpath('tests'); report=run_tests();
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    if nargin < 1, output_dir = fullfile(root,'results'); end
    if nargin < 2, make_plots = false; end
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    opts = struct('verbose',false,'residual_tolerance',1e-10);
    Nx = 13; Ny = 9; h = 2/(Nx-1); p = 1/(Ny-1);
    z = zeros(Nx,1);
    % Constant fields include zero Dirichlet values: no division by boundary data.
    for c = [0, 0.4]
        side = c*ones(Ny,1);
        initial = Initial_guess(0.2,h,p,Nx,Ny,z,z,side,side);
        [u,s] = Nonlinear_iteration(initial,0.2,h,p,Nx,Ny,z,z,side,side,500,1e-11,opts);
        assert(s.converged && norm(u-c,inf)<1e-10,'Constant solution regression.');
    end
    % An affine field exercises nonzero Neumann data and rectangular indexing.
    [X,Y] = ndgrid(linspace(0,2,Nx),linspace(0,1,Ny));
    exact = 0.2+0.1*X+0.05*Y;
    g = 0.05*ones(Nx,1); left = exact(1,:)'; right = exact(end,:)';
    opts.forcing = 0.1*exact;
    initial = Initial_guess(0.2,h,p,Nx,Ny,g,g,left,right);
    assert(norm(initial-exact(:),inf)<1e-11,'Initial-guess Neumann stencil regression.');
    % Nontrivial initial perturbation catches iterate layout errors after step one.
    initial = initial + 0.02*sin(pi*X(:)/2).*cos(pi*Y(:));
    [u,s] = Nonlinear_iteration(initial,0.2,h,p,Nx,Ny,g,g,left,right,500,1e-11,opts);
    assert(s.converged && s.iterations>1 && norm(u-exact(:),inf)<1e-9,'Affine field regression.');
    assert(s.dirichlet_error<1e-12 && s.neumann_error<1e-9,'Boundary regression.');
    % Evaluate residual independently using explicit stencils, not assembly code.
    assert(independent_residual(u,0.2,h,p,Nx,Ny,g,g,opts.forcing)<2e-10, ...
        'Independent discrete residual regression.');
    % A tiny relaxation must not create false convergence from a small update.
    saved_warning = warning('query','pde:NotConverged');
    cleanup = onCleanup(@() warning(saved_warning.state,'pde:NotConverged'));
    warning('off','pde:NotConverged');
    slow = opts; slow.relaxation = 1e-12;
    lastwarn('');
    [~,failed] = Nonlinear_iteration(initial,0.2,h,p,Nx,Ny,g,g,left,right,1,1e-6,slow);
    [~,warning_id] = lastwarn(); %#ok<ASGLU>
    assert(~failed.converged && failed.relative_update(end)<1e-6 && ...
        failed.residual(end)>slow.residual_tolerance,'False convergence regression.');
    clear cleanup;
    assert_error(@() Initial_guess(0,h,p,Nx,Ny,g,g,left,right));
    assert_error(@() Nonlinear_iteration(initial,0.2,h,p,Nx,Ny,g,g,left,right,10,1e-6,struct('relaxation',0)));
    assert_error(@() Initial_guess(0.2,h,p,Nx,Ny,g(1:end-1),g,left,right));
    % Fixed-epsilon refinement: nonlinear analytic 1-D solution and true 2-D MMS.
    analytic = refinement_study('analytic');
    manufactured = refinement_study('manufactured');
    assert(all(analytic.order(2:end)>1.8) && all(analytic.order(2:end)<2.2), ...
        'Analytic benchmark must converge at second order.');
    assert(all(manufactured.order(2:end)>1.8) && all(manufactured.order(2:end)<2.2), ...
        'Manufactured solution must converge at second order.');
    assert(manufactured.neumann_error(end)<manufactured.neumann_error(1)/10, ...
        'Independent Neumann derivative error must decrease under refinement.');
    report = struct('runtime',version,'analytic',analytic,'manufactured',manufactured);
    save(fullfile(output_dir,'validation.mat'),'report');
    write_csv(fullfile(output_dir,'analytic.csv'),analytic);
    write_csv(fullfile(output_dir,'manufactured.csv'),manufactured);
    if make_plots
        f = figure('Visible','off');
        loglog(analytic.h,analytic.error,'o-',manufactured.h,manufactured.error,'s-');
        xlabel('h_x'); ylabel('Maximum nodal error'); grid on;
        legend('1-D analytic solution','2-D manufactured solution','Location','northwest');
        print(f,fullfile(output_dir,'grid_convergence.png'),'-dpng','-r150'); close(f);
    end
    fprintf('All numerical regression tests passed. Runtime: %s\n',version);
end

function study = refinement_study(kind)
    nx = [21,41,81]; ny = [13,25,49];
    fields = {'h','error','order','residual','iterations','neumann_error'};
    for k = 1:numel(fields), study.(fields{k}) = zeros(3,1); end
    study.Nx = nx(:); study.Ny = ny(:); study.order(1) = NaN;
    for level = 1:3
        Nx = nx(level); Ny = ny(level); Lx = 2; Ly = 1;
        h = Lx/(Nx-1); p = Ly/(Ny-1); epsilon = 0.2;
        [X,Y] = ndgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny));
        if strcmp(kind,'analytic')
            % u=-a*tanh(a*(x-Lx/2)/(2*epsilon)) solves u*ux-epsilon*uxx=0.
            a = 0.5; exact = -a*tanh(a*(X-Lx/2)/(2*epsilon));
            forcing = zeros(Nx,Ny); bottom = zeros(Nx,1); top = bottom;
        else
            % Genuine y dependence, with unequal nonzero +y boundary derivatives.
            exact = 0.2+0.1*X+0.05*Y+0.02*Y.^2+0.03*sin(pi*X/Lx).*cos(pi*Y/Ly);
            ux = 0.1+0.03*(pi/Lx)*cos(pi*X/Lx).*cos(pi*Y/Ly);
            uxx = -0.03*(pi/Lx)^2*sin(pi*X/Lx).*cos(pi*Y/Ly);
            uyy = 0.04-0.03*(pi/Ly)^2*sin(pi*X/Lx).*cos(pi*Y/Ly);
            forcing = exact.*ux-epsilon*uxx-uyy;
            bottom = 0.05*ones(Nx,1); top = 0.09*ones(Nx,1);
        end
        left = exact(1,:)'; right = exact(end,:)';
        initial = Initial_guess(epsilon,h,p,Nx,Ny,bottom,top,left,right);
        opts = struct('forcing',forcing,'verbose',false,'residual_tolerance',1e-10);
        [u,info] = Nonlinear_iteration(initial,epsilon,h,p,Nx,Ny,bottom,top,left,right,1000,1e-12,opts);
        assert(info.converged,'Benchmark nonlinear solve failed.');
        independent = independent_residual(u,epsilon,h,p,Nx,Ny,bottom,top,forcing);
        assert(independent<2e-10,'Independent residual check failed.');
        assert(info.dirichlet_error<1e-12,'Dirichlet condition failed.');
        study.h(level) = h; study.error(level) = norm(u-exact(:),inf);
        study.residual(level) = independent; study.iterations(level) = info.iterations;
        study.neumann_error(level) = info.neumann_error;
        if level>1, study.order(level) = log2(study.error(level-1)/study.error(level)); end
        fprintf('%s %dx%d: error %.6e, order %.4f, residual %.3e\n', ...
            kind,Nx,Ny,study.error(level),study.order(level),independent);
    end
end

function value = independent_residual(u,epsilon,h,p,Nx,Ny,bottom,top,f)
    U = reshape(u,Nx,Ny); R = zeros(Nx-2,Ny);
    for j = 1:Ny
        center = U(2:end-1,j);
        ux = (U(3:end,j)-U(1:end-2,j))/(2*h);
        uxx = (U(3:end,j)-2*center+U(1:end-2,j))/h^2;
        if j==1
            ghost = U(2:end-1,2)-2*p*bottom(2:end-1);
            uyy = (U(2:end-1,2)-2*center+ghost)/p^2;
        elseif j==Ny
            ghost = U(2:end-1,end-1)+2*p*top(2:end-1);
            uyy = (ghost-2*center+U(2:end-1,end-1))/p^2;
        else
            uyy = (U(2:end-1,j+1)-2*center+U(2:end-1,j-1))/p^2;
        end
        R(:,j) = center.*ux-epsilon*uxx-uyy-f(2:end-1,j);
    end
    value = norm(R(:),inf);
end

function assert_error(callback)
    failed = false;
    try, callback(); catch, failed = true; end
    assert(failed,'Invalid input was accepted.');
end

function write_csv(path,study)
    fid = fopen(path,'w');
    if fid<0, error('pde:Output','Cannot write %s',path); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'Nx,Ny,h,max_error,observed_order,pde_residual,iterations,neumann_error\n');
    for k = 1:numel(study.h)
        fprintf(fid,'%d,%d,%.12g,%.12g,%.8g,%.12g,%d,%.12g\n', ...
            study.Nx(k),study.Ny(k),study.h(k),study.error(k),study.order(k), ...
            study.residual(k),study.iterations(k),study.neumann_error(k));
    end
end
