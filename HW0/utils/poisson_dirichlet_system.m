function [A, F_vec, X_full, Y_full, U_bc] = poisson_dirichlet_system(Nx, Ny, L, u_exact, f)
%POISSON_DIRICHLET_SYSTEM Build Dirichlet Poisson system on [0,L]^2.
%
% Constructs the 5-point Laplacian stencil for interior points and applies
% Dirichlet boundary conditions using the exact solution u_exact.
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
%   U_bc    - exact boundary values

h = L / (Nx + 1);   % uniform discretization
x_all  = linspace(0, L, Nx+2); % discretized x-grid
y_all  = linspace(0, L, Ny+2); % discretized y-grid
x = x_all(2:end-1); % grab interior points
y = y_all(2:end-1);
[X, Y] = meshgrid(x, y); % full 2D grid

% RHS evaluated on grid, scaled by h^2 for discrete Laplacian
F_grid = h^2 * f(X, Y);

% Dirichlet boundary data, corresponding to u_exact at boundaries
g_left   = @(y) sin(0).*cos(y); % x=0
g_right  = @(y) sin(L).*cos(y); % x=L
g_bottom = @(x) sin(x).*cos(0); % y=0
g_top    = @(x) sin(x).*cos(L); % y=L

% storage for spare matrix assembly
row_idx = []; col_idx = []; val_idx = [];

for j = 1:Ny
    for i = 1:Nx
        k = (j-1)*Nx + i;    % convert (i,j) grid index to linear index k
        current_x = X(j, i);
        current_y = Y(j, i);

        % center coefficient of Laplacian stencil
        row_idx = [row_idx; k]; col_idx = [col_idx; k]; val_idx = [val_idx; -4];

        % x-direction neighbors
        if i > 1    % left neighbor (i-1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
        else
            % Dirichlet BC at x=0 - set RHS
            F_grid(j,i) = F_grid(j,i) - g_left(current_y);
        end

        if i < Nx   % right neighbor (i+1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
        else
            % Dirichlet BC at x=L - set RHS
            F_grid(j,i) = F_grid(j,i) - g_right(current_y);
        end

        % y-direction neighbors
        if j > 1    % bottom neighbor (j-1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Nx]; val_idx = [val_idx; 1];
        else
            % Dirichlet BC at y=0 - set RHS
            F_grid(j,i) = F_grid(j,i) - g_bottom(current_x);
        end

        if j < Ny   % top neighbor (j+1)
            % interior neighbor
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Nx]; val_idx = [val_idx; 1];
        else
            % Dirichlet BC at y=L - set RHS
            F_grid(j,i) = F_grid(j,i) - g_top(current_x);
        end
    end
end

% build sparse matrix and RHS vector
A = sparse(row_idx, col_idx, val_idx, Nx*Ny, Nx*Ny);
F = F_grid.';
F_vec = F(:);

% put together grid for BCs
[X_full, Y_full] = meshgrid(x_all, y_all);
U_bc = u_exact(X_full, Y_full);
end


