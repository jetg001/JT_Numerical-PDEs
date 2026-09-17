%% Problem 1 - First Derivative Approximation
% This script computes the first derivative of a smooth function using
% centered finite differences. Convergence is measured in the l_inf and l_2
% norms, and the numerical solution is compared against the exact
% derivative.
%
% Utilities used:
%   - first_derivative_fd.m
%   - format_loglog_axes.m
%
% Figures generated:
%   - p1_derivative.png
%   - p1_convergence.png

clear; close all; clc;
addpath('utils\');

f  = @(x) exp(sin(x));
fp = @(x) cos(x).*exp(sin(x));

N  = 200;
x  = linspace(0, 2*pi, N);
dx = x(2) - x(1);
fx = f(x);

fp_num   = first_derivative_fd(fx, dx);
fp_exact = fp(x);

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(x, fp_exact, 'Color', [0 0.6275 0.8431], 'LineWidth', 10); hold on;
plot(x, fp_num, '--', 'Color', [1 0.5 0], 'LineWidth', 14)

xlim([0 2*pi]);
ylim([min(fp_num)-0.1, max(fp_num)+0.1]);

set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi], 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [-1.5, -1, -0.5, 0, 0.5, 1, 1.5], 'YTickLabel', {'$-1.5$', '$-1$', '$-0.5$', '$0$', '$0.5$', '$1$', '$1.5$'}, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 30;

xlabel('$x$', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Derivative', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Exact vs Numerical Derivative of $f(x)=\exp(\sin(x))$', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend({'Exact $f''(x)$', 'Numerical $f''_n(x)$'},'FontSize', 30, 'FontWeight', 'bold', 'Location', 'north', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p1_derivative.png'));
hold off;


Ns = round(logspace(2, 4, 20));
err_inf = zeros(size(Ns));
err_2   = zeros(size(Ns));

for k = 1:length(Ns)
    Nk = Ns(k);
    xk = linspace(0, 2*pi, Nk);
    dxk = xk(2) - xk(1);
    fxk = f(xk);
    fp_num_k   = first_derivative_fd(fxk, dxk);
    fp_exact_k = fp(xk);

    err_inf(k) = norm(fp_exact_k - fp_num_k, Inf) / norm(fp_exact_k, Inf);
    err_2(k)   = norm(fp_exact_k - fp_num_k, 2)   / norm(fp_exact_k, 2);
end

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
loglog(Ns, err_inf, 'o-', 'Color', [0 0.6275 0.8431], 'LineWidth', 10, 'MarkerSize', 10);
hold on;
loglog(Ns, err_2, 's-', 'Color', [0.81 0.47 0.66], 'LineWidth', 10, 'MarkerSize', 10);
loglog(Ns, Ns.^(-1), '--', 'Color', [0.83 0.37 0], 'LineWidth', 10);
loglog(Ns, Ns.^(-3/2), '--', 'Color', [1 0.5 0],'LineWidth', 10);
format_loglog_axes();
xlabel('$N$ gridpoints', 'Interpreter', 'latex');
ylabel('Relative Error', 'Interpreter', 'latex');
title('$\ell_\infty$ vs $\ell_2$ Error for Numerical Derivative of $f(x)=\exp(\sin(x))$', 'Interpreter', 'latex');
legend({'$\ell_\infty$ Error', '$\ell_2$ Error', 'Reference Slope $N^{-1}$', 'Reference Slope $N^{-3/2}$'}, 'Location', 'southwest', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p1_convergence.png'));
hold off;


