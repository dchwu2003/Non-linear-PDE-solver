%% The following code snippet uses finite difference method and
% interative method to approximate the solution of 2-D second-order 
% nonlinear PDE U*Ux-Uyy-εUxx=0 with the boundary
%conditions given. 

clear; clc; close all;

%Set up a sequence of epsilon
N_epsilon = 5;
Epsilon = linspace(0.1,0.05,N_epsilon);
% Total number of nodes in the 2-D grid
Nx = 100; Ny = 100;  
%set total range of spanning of the grid
spanx = 6; spany = 5;
%Set up the spacing 
h = spanx/(Nx-1); 
p = spany/(Ny-1);

% Check if the grid is dense enough
if h > (2*Epsilon(end)) || p > (2*Epsilon(end))
   disp('The grid is not dense enough. Program stops.');
   return;
end    

%set the parameters for the iteration process
max_iter = 30;
err_ctrl = 1e-3;
prev_sol = [];
stablize_early = false;
U_iter = [];

%Input the boundary conditions
inibotbd = linspace(0,spanx, Nx);
initopbd = linspace(0,spanx, Nx);
inileftbd = linspace(0,spany, Ny);
inirightbd = linspace(0,spany, Ny);

leftbd = ones(length(inileftbd),1);
rightbd = -1.1*ones(length(inirightbd),1);
upperbd = zeros(length(initopbd),1);
bottombd = zeros(length(inibotbd),1);

 % Step1, solve εUxx+Uyy=0 to obtain the inital guess
U_initial = Initial_guess(Epsilon(1), h, p, Nx, Ny, ...
                              bottombd, upperbd,leftbd,rightbd);

% Step2, use the inital guess to start the iterative processive of
% V*Ux-εUxx-Uyy=0 where V is given by the previous iteration.
for iterator = 1:length(Epsilon)
  
    if iterator == 1
        U_iter = U_initial;

    else 
        U_iter = curr_sol;
    end
 

    curr_sol = Nonlinear_iteration(U_iter, Epsilon(iterator), h, ...
                       p, Nx, Ny, bottombd, upperbd,leftbd, ...
                                  rightbd, max_iter, err_ctrl);
    
    %Graph the solution per X times iteration
    if mod(iterator, 2) == 1 && iterator ~= length(Epsilon)  
        Mat_epsilon = (reshape(curr_sol,Nx,Ny))';
        figure;
        surf(linspace(0, spanx, Nx), linspace(0, spany, Ny),Mat_epsilon);
        title(sprintf('Graph of solution upon %dth iteration', iterator));
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('u(x,y)');
        shading interp
        colorbar;
    end

    %Convergence test for sequence of final solutions corresponding to each
    %epsilon
    if ~isempty(prev_sol)
     
        rel_err = norm(curr_sol - prev_sol) / norm(prev_sol);

        fprintf(['Relative error from between two solutions with respect to' ...
            ' two different epsilon is %.4e\n'], rel_err);
       
        if rel_err < err_ctrl && iterator < N_epsilon
            fprintf('Solution stabilizing early\n');
            stablize_early = true;
            break;
        end
    else 
        prev_sol = curr_sol;
    end


end

if stablize_early
    fprintf('Converged early at %dth epsilon', iterator);
else
    fprintf('Completed all epsilon values\n');
end

%graph the initial guess
figure;
Mat_Initial = (reshape(U_initial,Nx,Ny))';
surf(linspace(0, spanx, Nx), linspace(0, spany, Ny), Mat_Initial);
title('Initial Guess');
xlabel('x-axis');
ylabel('y-axis');
zlabel('u(x,y)');
shading interp
colorbar;


%graph the ultimate result
figure;
Mat_final = (reshape(curr_sol,Nx,Ny))'; 
surf(linspace(0, spanx, Nx), linspace(0, spany, Ny), Mat_final);
title('Approximate Solution of UUx-U_{yy} = 0 in 2-D');
xlabel('x-axis');
ylabel('y-axis');
zlabel('u(x,y)');
shading interp
colorbar;