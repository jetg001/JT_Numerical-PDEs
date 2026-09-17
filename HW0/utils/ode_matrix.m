function [A, f_i] = ode_matrix(x_i, dx, u0, uN, f)
%ODE_MATRIX Build sparse matrix for u'' + sin(x)u' + u = f(x).
%
% Constructs the tridiagonal matrix corresponding to the centered finite
% differences on interior points x_i. Dirichlet boundary conditions are
% applied at endpoints.
%
% Inputs:
%   x_i - interior grid points
%   dx  - grid spacing
%   u0  - left boundary value
%   uN  - right boundary value
%   f   - forcing function
%
% Outputs:
%   A   - sparse tridiagonal matrix
%   f_i - right-hand side vector

N = numel(x_i);
s_i = sin(x_i);

a_i = 1/dx^2 - s_i/(2*dx);
b_i = -2/dx^2 + 1;
c_i = 1/dx^2 + s_i/(2*dx);

rows = [2:N, 1:N, 1:N-1];
cols = [1:N-1, 1:N, 2:N];
vals = [a_i(2:end), b_i*ones(1,N), c_i(1:end-1)];

A = sparse(rows, cols, vals, N, N);

f_i      = f(x_i)';
f_i(1)   = f_i(1)   - a_i(1)*u0;
f_i(end) = f_i(end) - c_i(end)*uN;
end