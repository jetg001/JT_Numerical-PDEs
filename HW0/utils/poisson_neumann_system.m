function [A, F_vec, X, Y] = poisson_neumann_system(Nx, Ny, L, u_exact, f)
%POISSON_NEUMANN_SYSTEM Build Neumann Poisson system on [0,L]^2.
%
% Constructs the 5-point Laplacian stencil for interior points and enforces
% Neumann boundary conditions using ghost points. The resulting matrix has
% a one-dimensional nullspace and must be pinned externally.
%
% Inputs:
%   Nx, Ny  - number of interior points in x and y
%   L       - domain length
%   u_exact - exact solution function
%   f       - forcing term (RHS)
%
% Outputs:
%   A       - sparse matrix
%   F_vec   - right-hand side vector
%   X_full  - full grid (including boundaries)
%   Y_full  - full grid (including boundaries)

h = L / (Nx - 1);       % uniform discretization
x = linspace(0, L, Nx); % discretized x-grid
y = linspace(0, L, Ny); % discretized y-grid
[X,Y] = meshgrid(x, y); % full 2D grid

% RHS evaluated on grid, scaled by h^2 for discrete Laplacian
F_grid = h^2 * f(X,Y);

% Neumann boundary data, corresponding to normal vectors at boundaries
g_left   = @(y) -cos(0)*cos(y); % x = 0
g_right  = @(y) cos(L)*cos(y);  % x = L
g_bottom = @(x) 0*x;            % y = 0
g_top    = @(x) -sin(x)*sin(L); % y = L

% storage for sparse matrix assembly
row_idx = []; col_idx = []; val_idx = [];

for j = 1:Ny
    for i = 1:Nx
        k = (i-1)*Ny + j;   % convert (i,j) grid index to linear index k
        current_x = X(j,i);
        current_y = Y(j,i);

        % center coefficient of Laplacian stencil
        row_idx = [row_idx; k]; col_idx = [col_idx; k]; val_idx = [val_idx; -4];

        % x-direction neighbors
        if i > 1    % left neighbor (i-1)
            % interior neigbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Ny]; val_idx = [val_idx; 1];
        else
            % Neumann BC at x=0: use ghost point and mirror right neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Ny]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_left(current_y); % modify RHS
        end

        if i < Nx   % right neighbor (i+1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Ny]; val_idx = [val_idx; 1];
        else
            % Neumann BC at x=L: use ghost point and mirror left neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Ny]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_right(current_y);
        end

        % y-direction neighbors
        if j > 1    % bottom neighbor (j-1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
        else
             % Neumann BC at y=0: use ghost point and mirror top neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_bottom(current_x);
        end

        if j < Ny   % top neighbor (j+1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
        else
            % Neumann BC at y=L: use ghost point and mirror bottom neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_top(current_x);
        end
    end
end

% build sparse matrix and RHS vector
A = sparse(row_idx, col_idx, val_idx, Nx*Ny, Nx*Ny);
F_vec = F_grid(:);

% Neumann Laplacian singular -> pin one degree of freedom
pin = 1;         % choose first grid point
A(pin, :) = 0;   % zero out row
A(pin, pin) = 1; % enforce u(pin) = u_exact
F_vec(pin) = u_exact(X(pin), Y(pin));
end



