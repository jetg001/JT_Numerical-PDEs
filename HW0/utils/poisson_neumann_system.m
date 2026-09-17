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

h = L / (Nx - 1);
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
[X,Y] = meshgrid(x, y);

F_grid = h^2 * f(X,Y);

g_left   = @(y) -cos(0)*cos(y);
g_right  = @(y) cos(L)*cos(y);
g_bottom = @(x) 0*x;
g_top    = @(x) -sin(x)*sin(L);

row_idx = []; col_idx = []; val_idx = [];

for j = 1:Ny
    for i = 1:Nx
        k = (i-1)*Ny + j;
        current_x = X(j,i);
        current_y = Y(j,i);

        row_idx = [row_idx; k]; col_idx = [col_idx; k]; val_idx = [val_idx; -4];

        if i > 1
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Ny]; val_idx = [val_idx; 1];
        else
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Ny]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_left(current_y);
        end

        if i < Nx
            row_idx = [row_idx; k]; col_idx = [col_idx; k+Ny]; val_idx = [val_idx; 1];
        else
            row_idx = [row_idx; k]; col_idx = [col_idx; k-Ny]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_right(current_y);
        end

        if j > 1
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
        else
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_bottom(current_x);
        end

        if j < Ny
            row_idx = [row_idx; k]; col_idx = [col_idx; k+1]; val_idx = [val_idx; 1];
        else
            row_idx = [row_idx; k]; col_idx = [col_idx; k-1]; val_idx = [val_idx; 1];
            F_grid(j,i) = F_grid(j,i) - 2*h*g_top(current_x);
        end
    end
end

A = sparse(row_idx, col_idx, val_idx, Nx*Ny, Nx*Ny);
F_vec = F_grid(:);

pin = 1;
A(pin, :) = 0;
A(pin, pin) = 1;
F_vec(pin) = u_exact(X(pin), Y(pin));
end



