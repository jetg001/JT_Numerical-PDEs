%% Problem 3 - ODE: u'' + sin(x)u' + u = f(x)
% This script solves a linear ODE with variable coefficients using centered
% finite differences. Dirichlet boundary conditinos are applied at x = 0
% and x = 5. Convergence is measured in the l_inf and l_2 norms.
%
% Utilities used:
%   - ode_matrix.m
%   - format_loglog_axes.m
%
% Figures generated:
%   - p3_solution.png
%   - p3_convergence.png

clear; close all; clc;
addpath('utils\');

u_exact = @(x) sin(x);
f       = @(x) sin(x).*cos(x);

a = 0; b = 5;
u0 = u_exact(a);
uN = u_exact(b);

Ns = round(logspace(2, 4, 20));
err_inf = zeros(size(Ns));
err_2   = zeros(size(Ns));

for k = 1:length(Ns)
    N  = Ns(k);
    x  = linspace(a, b, N+1);
    dx = x(2) - x(1);

    x_i = x(2:N);
    [A, f_i] = ode_matrix(x_i, dx, u0, uN, f);

    u_sol = A \ f_i;
    u_num = [u0; u_sol; uN];
    u_ex  = u_exact(x)';

    err_inf(k) = norm(u_ex - u_num, Inf) / norm(u_ex, Inf);
    err_2(k)   = norm(u_ex - u_num, 2)   / norm(u_ex, 2);

    if k == 1
        figure;
        plot(x, u_ex, 'Color', [0 0.6275 0.8431], 'LineWidth', 18); hold on;
        plot(x, u_num, '--', 'Color', [1 0.5 0], 'LineWidth', 25);
        ylim([min(u_num)-0.1, max(u_num)+0.1]);
        set(gca, 'XTick', 0:5, 'XTickLabel', {'$0$', '$1$', '$2$', '$3$', '$4$', '$5$'}, 'TickLabelInterpreter', 'latex');
        ax = gca; ax.FontSize = 35;
        xlabel('$x$', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
        ylabel('Solution to ODE', 'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
        title('Exact vs Numerical Solution to $u^{\prime\prime}+\sin(x)u^\prime+u(x)=f(x)$','FontSize', 35, 'FontWeight', 'bold', ...
            'Interpreter', 'latex');
        legend({'Exact: $u(x) = \sin(x)$', 'Numerical: $u_n(x)$'}, 'FontSize', 45, 'FontWeight', 'bold', 'Location', 'southwest', 'Interpreter', 'latex');
        saveas(gcf, fullfile('figures', 'p3_solution.png'));
        hold off;
    end
end

figure;

loglog(Ns, err_inf, 'o-', 'Color', [0 0.6275 0.8431], 'LineWidth', 11, 'MarkerSize', 13); hold on;
loglog(Ns, err_2, 's-', 'Color', [0.81 0.47 0.66], 'LineWidth', 9, 'MarkerSize', 11);
grid on;
loglog(Ns, Ns.^(-2), '--', 'Color', [0.83 0.37 0], 'LineWidth', 10);
set(gca, 'XTick', [1e2 1e3 1e4], 'XTickLabel', {'$10^2$', '$10^3$', '$10^4$'}, 'TickLabelInterpreter', 'latex');
format_loglog_axes();
xlabel('$N$ gridpoints','Interpreter', 'latex');
ylabel('Relative Error', 'Interpreter', 'latex');
title('$\ell_\infty$ vs $\ell_2$ Error for Numerical ODE Solution', 'Interpreter', 'latex');
legend({'$\ell_\infty$  Error', '$\ell_2$ Error', 'Reference Slope $N^{-2}$'}, 'Location', 'northeast', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p3_convergence.png'));
hold off;
