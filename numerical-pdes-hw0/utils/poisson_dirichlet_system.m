function [A, F_vec, X_full, Y_full, U_bc] = poisson_dirichlet_system(Nx, Ny, L, u_exact)
%POISSON_DIRICHLET_SYSTEM Build Dirichlet Poisson system on [0,L]^2.
%
% Constructs the 5-point Laplacian stencil for interior points and applies
% Dirichlet boundary conditions using the exact solution u_exact.
%
% Inputs:
%   Nx, Ny  - number of interior points in x and y
%   L       - domain length
%   u_exact - exact solution function
%
% Outputs:
%   A       - sparse matrix
%   F_vec   - right-hand side vector
%   X_full  - full grid (including boundaries)
%   Y_full  - full grid (including boundaries)
%   U_bc    - exact boundary values

h = L / (Nx + 1);
x_all  = linspace(0, L, Nx+2);
y_all  = linspace(0, L, Ny+2);
x = x_all(2:end-1);
y = y_all(2:end-1);
[X, Y] = meshgrid(x, y);

f = @(x,y) -2*sin(x).*cos(y);
F_grid = h^2 * f(X, Y);

g_left   = @(y) sin(0).*cos(y);
g_right  = @(y) sin(L).*cos(y);
g_bottom = @(x) sin(x).*cos(0);
g_top    = @(x) sin(x).*cos(L);

row_idx = []; col_idx = []; val_idx = [];
for j = 1:Ny
    for i = 1:Nx
        k = (j-1)*Nx + i;
        current_x = X(j, i);
        current_y = Y(j, i);

        row_idx = [row_idx; k]; col_idx = [col_idx; k]; val_idx = [val_idx; -4];

        if i > 1
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
        else
            F_grid(j,i) = F_grid(j,i) - g_left(current_y);
        end

        if i < Nx
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
        else
            F_grid(j,i) = F_grid(j,i) - g_right(current_y);
        end

        if j > 1
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Nx]; val_idx = [val_idx; 1];
        else
            F_grid(j,i) = F_grid(j,i) - g_bottom(current_x);
        end

        if j < Ny
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Nx]; val_idx = [val_idx; 1];
        else
            F_grid(j,i) = F_grid(j,i) - g_top(current_x);
        end
    end
end

A = sparse(row_idx, col_idx, val_idx, Nx*Ny, Nx*Ny);
F = F_grid.';
F_vec = F(:);

[X_full, Y_full] = meshgrid(x_all, y_all);
U_bc = u_exact(X_full, Y_full);
end


