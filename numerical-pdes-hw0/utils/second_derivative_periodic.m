function f2p_num = second_derivative_periodic(fx, dx)
%SECOND_DERIVATIVE_PERIODIC    Periodic centered finite difference second
% derivative.
%
% f2p_num = SECOND_DERIVATIVE_PERIODIC(fx, dx) computes the second
% derivative using the centered finite difference formula with periodic
% wrapping:
%       f''(x_i) ≈ (f_{i+1} -2 f_i + f_{i-1}) / dx^2 
% Inputs: 
%   fx - function values at grid points 
%   dx - grid spacing 
% Output: 
%   f2p_num - numerical first derivative

N = numel(fx);
f2p_num = zeros(size(fx));

f2p_num(1) = (fx(2) - 2*fx(1) + fx(N-1)) / dx^2;
f2p_num(N) = f2p_num(1);
f2p_num(2:N-1) = (fx(3:N) - 2*fx(2:N-1) + fx(1:N-2)) / dx^2;
end