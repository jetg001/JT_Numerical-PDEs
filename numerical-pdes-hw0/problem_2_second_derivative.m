%% Problem 2 - Second Derivative with Periodic Boundary Conditions
% This script computes the second derivative of a periodic function using
% centered finite differences with periodic wrapping. Convergence is 
% measured in the l_inf and l_2 norms, and the numerical solution is 
% compared against the exact second derivative.
%
% Utilities used:
%   - second_derivative_fd.m
%   - format_loglog_axes.m
%
% Figures generated:
%   - p2_second_derivative.png
%   - p2_convergence.png

clear; close all; clc;
addpath('utils\');

f   = @(x) exp(sin(x));
f2p = @(x) exp(sin(x)).*(cos(x).^2 - sin(x));

N = 200;
x = linspace(0, 2*pi, N);
dx = x(2) - x(1);
fx = f(x);

f2p_num   = second_derivative_periodic(fx, dx);
f2p_exact = f2p(x);

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(x, f2p_exact, 'Color', [0 0.6275 0.8431], 'LineWidth', 10);
hold on;
plot(x, f2p_num, '--', 'Color', [1 0.5 0], 'LineWidth', 14);

xlim([0 2*pi]); ylim([min(f2p_num)-0.1, max(f2p_num)+0.1]);
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi], 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = 30;

xlabel('$x$', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Second Derivative', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Exact vs Numerical Second Derivative of $f(x)=\exp(\sin(x))$','FontSize', 30, 'FontWeight', 'bold', ...
    'Interpreter', 'latex');
legend({'Exact: $f^{\prime\prime}(x)$', 'Numerical: $f_n^{\prime\prime}(x)$'}, 'FontSize', 30, 'FontWeight', 'bold', 'Location', 'east', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p2_second_derivative.png'));
hold off;


Ns = round(logspace(2, 4, 20));
err_inf = zeros(size(Ns));
err_2   = zeros(size(Ns));

for k = 1:length(Ns)
    Nk = Ns(k);
    xk = linspace(0, 2*pi, Nk);
    dxk = xk(2) - xk(1);
    fxk = f(xk);
    f2p_num_k   = second_derivative_periodic(fxk, dxk);
    f2p_exact_k = f2p(xk);

    err_inf(k) = norm(f2p_exact_k - f2p_num_k, Inf) / norm(f2p_exact_k, Inf);
    err_2(k)   = norm(f2p_exact_k - f2p_num_k, 2)   / norm(f2p_exact_k, 2);
end

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
loglog(Ns, err_inf, 'o-', 'Color', [0 0.6275 0.8431], 'LineWidth', 12, 'MarkerSize', 12);
hold on;
loglog(Ns, err_2, 's-', 'Color', [0.81 0.47 0.66], 'LineWidth', 10, 'MarkerSize', 10);
loglog(Ns, Ns.^(-2), '--', 'Color', [0.83 0.37 0], 'LineWidth', 10);
format_loglog_axes();
xlabel('$N$ gridpoints', 'Interpreter', 'latex');
ylabel('Relative Error', 'Interpreter', 'latex');
title('$\ell_\infty$ vs $\ell_2$ Error for Numerical Derivative of $f(x)=\exp(\sin(x))$', 'Interpreter', 'latex');
legend({'$\ell_\infty$ Error', '$\ell_2$ Error', 'Reference Slope $N^{-2}$'}, 'Location', 'southwest', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p2_convergence.png'));
hold off;