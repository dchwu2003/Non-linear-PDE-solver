function result = run_example(output_dir, make_plots)
%RUN_EXAMPLE Fixed-grid epsilon continuation, saving parameters and diagnostics.
% result=run_example() creates results/ relative to this file.
% run_example(tempdir,false) runs without graphics.
    if nargin < 1, output_dir = fullfile(fileparts(mfilename('fullpath')),'results'); end
    if nargin < 2, make_plots = true; end
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    config = default_config();
    x = linspace(0,config.Lx,config.Nx); y = linspace(0,config.Ly,config.Ny);
    h = x(2)-x(1); p = y(2)-y(1);
    left = config.left*ones(config.Ny,1); right = config.right*ones(config.Ny,1);
    bottom = zeros(config.Nx,1); top = bottom;
    u = Initial_guess(config.epsilon(1),h,p,config.Nx,config.Ny,bottom,top,left,right);
    initial = u;
    solutions = zeros(config.Nx*config.Ny,numel(config.epsilon));
    history = cell(numel(config.epsilon),1);
    epsilon_change = NaN(numel(config.epsilon),1);
    options = struct('relaxation',config.relaxation,'verbose',false, ...
        'residual_tolerance',config.residual_tolerance);
    for k = 1:numel(config.epsilon)
        previous = u;
        [u, info] = Nonlinear_iteration(u,config.epsilon(k),h,p,config.Nx,config.Ny, ...
            bottom,top,left,right,config.max_iterations,config.update_tolerance,options);
        if ~info.converged, error('pde:ExampleNotConverged','Failed at epsilon=%g.',config.epsilon(k)); end
        if k > 1, epsilon_change(k) = norm(u-previous,inf)/max(1,norm(previous,inf)); end
        % Every requested epsilon is solved. A small parameter change is not a
        % convergence test for the PDE or a reason to skip the target epsilon.
        solutions(:,k) = u; history{k} = info;
        fprintf('epsilon=%.4f: %d iterations, residual %.3e, cell Pe %.3f\n', ...
            config.epsilon(k),info.iterations,info.residual(end),info.max_cell_peclet);
    end
    result = struct('config',config,'x',x,'y',y,'initial',initial, ...
        'solutions',solutions,'history',{history},'epsilon_change',epsilon_change);
    save(fullfile(output_dir,'example.mat'),'result');
    fid = fopen(fullfile(output_dir,'continuation.csv'),'w');
    if fid < 0, error('pde:Output','Cannot write continuation.csv'); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'epsilon,iterations,pde_residual,relative_epsilon_change,max_cell_peclet\n');
    for k = 1:numel(history)
        s = history{k}; fprintf(fid,'%.8g,%d,%.12g,%.12g,%.8g\n', ...
            config.epsilon(k),s.iterations,s.residual(end),epsilon_change(k),s.max_cell_peclet);
    end
    if make_plots
        f = figure('Visible','off');
        surf(x,y,reshape(u,config.Nx,config.Ny)'); shading interp; colorbar;
        xlabel('x'); ylabel('y'); zlabel('u'); title(sprintf('Final solution, epsilon = %.3f',config.epsilon(end)));
        print(f,fullfile(output_dir,'solution.png'),'-dpng','-r150'); close(f);
        f = figure('Visible','off'); hold on;
        for k = 1:numel(history), semilogy(history{k}.residual); end
        set(gca,'YScale','log'); grid on; xlabel('Picard iteration'); ylabel('PDE residual (infinity norm)');
        legend(arrayfun(@(e) sprintf('epsilon = %.4f',e),config.epsilon,'UniformOutput',false));
        print(f,fullfile(output_dir,'residuals.png'),'-dpng','-r150'); close(f);
        f = figure('Visible','off'); hold on;
        for k = 1:numel(history)
            M = reshape(solutions(:,k),config.Nx,config.Ny); plot(x,M(:,ceil(config.Ny/2)));
        end
        xlabel('x'); ylabel('u at mid-height'); grid on;
        legend(arrayfun(@(e) sprintf('epsilon = %.4f',e),config.epsilon,'UniformOutput',false));
        print(f,fullfile(output_dir,'profiles.png'),'-dpng','-r150'); close(f);
    end
end
