function [A, b, Dx, free] = assemble_system(epsilon, h, p, Nx, Ny, bottom, top, left, right, forcing)
%ASSEMBLE_SYSTEM Diffusion and centered x derivative on an x-fast grid.
% U is Nx-by-Ny and u=U(:); bottom/top prescribe +y derivatives, not normals.
% Neumann ghost values are eliminated using centered derivatives. Dirichlet
% conditions take precedence at the four physical corners.
    validateattributes(epsilon, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(h, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(p, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(Nx, {'numeric'}, {'scalar','integer','>=',3});
    validateattributes(Ny, {'numeric'}, {'scalar','integer','>=',3});
    bottom = boundary_vector(bottom, Nx, 'bottom');
    top = boundary_vector(top, Nx, 'top');
    left = boundary_vector(left, Ny, 'left');
    right = boundary_vector(right, Ny, 'right');
    if nargin < 10 || isempty(forcing), forcing = zeros(Nx, Ny); end
    validateattributes(forcing, {'numeric'}, {'real','finite','size',[Nx, Ny]});
    n = Nx*Ny;
    ids = reshape(1:n, Nx, Ny);
    free = ids(2:end-1, :); free = free(:);
    west = free-1; east = free+1;
    rx = epsilon/h^2; ry = 1/p^2;
    rows = [free; free; free];
    cols = [free; west; east];
    vals = [(2*rx+2*ry)*ones(size(free)); -rx*ones(size(free)); -rx*ones(size(free))];
    for j = 1:Ny
        r = ids(2:end-1, j);
        if j == 1
            rows = [rows; r]; cols = [cols; r+Nx]; vals = [vals; -2*ry*ones(size(r))]; %#ok<AGROW>
        elseif j == Ny
            rows = [rows; r]; cols = [cols; r-Nx]; vals = [vals; -2*ry*ones(size(r))]; %#ok<AGROW>
        else
            rows = [rows; r; r]; cols = [cols; r-Nx; r+Nx]; %#ok<AGROW>
            vals = [vals; -ry*ones(2*numel(r),1)]; %#ok<AGROW>
        end
    end
    bd = [ids(1,:)'; ids(end,:)'];
    A = sparse([rows; bd], [cols; bd], [vals; ones(size(bd))], n, n);
    Dx = sparse([free; free], [west; east], ...
        [-ones(size(free)); ones(size(free))]/(2*h), n, n);
    b = forcing(:);
    r = ids(2:end-1,1); b(r) = b(r) - 2*bottom(2:end-1)/p;
    r = ids(2:end-1,end); b(r) = b(r) + 2*top(2:end-1)/p;
    b(ids(1,:)) = left; b(ids(end,:)) = right;
end

function v = boundary_vector(v, count, name)
    validateattributes(v, {'numeric'}, {'vector','numel',count,'real','finite'}, mfilename, name);
    v = v(:);
end
