function fp_num = first_derivative_fd(fx, dx)
%FIRST_DERIVATIVE_FD   First derivative using centered + one-sided differences.
%
% fp = FIRST_DERIVATIVE_FD(fx, dx) computes the first derivative of a
% function samples at uniform grid points using the centered finite
% difference formula: 
%       f'(x_i) ≈ (f_{i+1} - f_{i-1}) / (2 dx) 
% Inputs: 
%   fx - function values at grid points 
%   dx - grid spacing 
% Output: 
%   fp - numerical first derivative

N = numel(fx);
fp_num = zeros(size(fx));

fp_num(1) = (fx(2) - fx(1)) / dx;
fp_num(N) = (fx(N) - fx(N-1)) / dx;
fp_num(2:N-1) = (fx(3:N) - fx(1:N-2)) / (2*dx);
end