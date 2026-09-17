%% Problem 4 - Poisson Equation with Dirichlet Boundary Conditions
% This script solves the 2D Poisson equation on [0,L'^2 using a 5-point
% Laplacian stencil with Dirichlet boundary conditions. The numerical
% solution is compared against the exact solution u(x,y) = sin(x)cos(y).
% Convergence is measured in the l_inf and l_2 norms.
%
% Utilities used:
%   - poisson_dirichlet_system.m
%   - slanCM.m
%   - format_loglog_axes.m
%
% Figures generated:
%   - p4_dirichlet_solution.png
%   - p4_dirichlet_convergence.png

clear; close all; clc;
addpath('utils\');

u_exact = @(x,y) sin(x).*cos(y);
L = 5;
Ns = round(logspace(1.5, 2.5, 12));
err_inf = zeros(size(Ns));
err_2   = zeros(size(Ns));

for m = 1:length(Ns)
    Nx = Ns(m); Ny = Nx;
    [A, F_vec, X_full, Y_full, U_bc] = poisson_dirichlet_system(Nx, Ny, L, u_exact);

    U_vec  = A \ F_vec;
    U_grid = reshape(U_vec, Ny, Nx).';

    U = U_bc;
    U(2:end-1, 2:end-1) = U_grid;

    u_ex = u_exact(X_full, Y_full);

    err_inf(m) = norm(u_ex - U, Inf) / norm(u_ex, Inf);
    err_2(m)   = norm(u_ex - U, 2)   / norm(u_ex, 2);

    if m == 1
        figure;
        colormap(slanCM('viridis'));
        cmin = min([U(:); u_ex(:)]);
        cmax = max([U(:); u_ex(:)]);
        t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
        ax1 = nexttile;
        surf(X_full,Y_full, U, 'EdgeColor','none'); view(2); shading interp; hold on;
        xlim([0 5]); ylim([0,5]); clim([cmin cmax]);
        set(gca, 'XTick', 0:5, 'YTick', 0:5, 'TickLabelInterpreter', 'latex')
        ax = gca; ax.FontSize = 25;
        xlabel('$x$', 'FontSize', 25, 'FontWeight', 'bold', 'Interpreter', 'latex');
        ylabel('$y$', 'FontSize', 25, 'FontWeight', 'bold', 'Interpreter', 'latex');
        title('Numerical: $u_n(x,y)$', 'FontSize', 25, 'FontWeight', 'bold', 'Interpreter', 'latex');
        
        ax2 = nexttile;
        surf(X_full,Y_full, u_ex, 'EdgeColor','none'); view(2); shading interp; hold on;
        xlim([0 5]); ylim([0 5]); clim([cmin cmax]);
        set(gca, 'XTick', [0 1 2 3 4 5], 'YTick', [0 1 2 3 4 5], 'TickLabelInterpreter', 'latex');
        ax = gca; ax.YAxisLocation = 'right'; ax.FontSize = 25;
        xlabel('$x$', 'FontSize', 25, 'FontWeight', 'bold', 'Interpreter', 'latex');
        title('Exact: $u(x,y)=\sin(x)\cos(y)$', 'FontSize', 25, 'FontWeight', 'bold', 'Interpreter', 'latex');
        title(t, 'Numerical (Dirichlet) vs Exact Solution of the Poisson Equation', 'FontSize', 28, 'FontWeight', 'bold', 'Interpreter', 'latex');
        
        cb = colorbar; cb.Parent = t; cb.FontSize = 20;  cb.TickLabelInterpreter = 'latex';
        drawnow;
        pos1 = ax1.Position; pos2 = ax2.Position;
        gap_start = pos1(1) + pos1(3); gap_end = pos2(1);
        cb_width = 0.025; cb_left = gap_start + (gap_end - gap_start)/2 - cb_width/2 + 0.05;
        height_factor = 1.35; new_height = pos1(4)*height_factor; 
        new_bottom = pos1(2) - (new_height - pos1(4))/2;
        cb.Position = [cb_left, new_bottom, cb_width, new_height+0.01];
        saveas(gcf, fullfile('figures', 'p4_dirichlet_solution.png'));
        hold off;
    end
end

figure;
loglog(Ns, err_inf, 'o-', 'Color', [0 0.6275 0.8431], 'LineWidth', 11, 'MarkerSize', 13); hold on;
loglog(Ns, err_2, 's-', 'Color', [0.81 0.47 0.66], 'LineWidth', 9, 'MarkerSize', 11);
grid on;
loglog(Ns, Ns.^(-2), '--', 'Color', [0.83 0.37 0], 'LineWidth', 10);
set(gca, 'XTick', [50 1e2 2e2 3e2], 'XTickLabel', { '$50$', '$100$', '$200$', '$300$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [1e-5 1e-4 1e-3], 'YTickLabel', { '$10^{-5}$', '$10^{-4}$', '$10^{-3}$'})
format_loglog_axes();
xlabel('$N$ gridpoints', 'Interpreter', 'latex');
ylabel('Relative Error', 'Interpreter', 'latex');
title('$\ell_\infty$ vs $\ell_2$ Error for Numerical Poisson (Dirichlet)', 'Interpreter', 'latex');
legend({'$\ell_\infty$ Error', '$\ell_2$ Error', 'Reference Slope $N^{-2}$'}, 'Location', 'northeast', 'Interpreter', 'latex');
saveas(gcf, fullfile('figures', 'p4_dirichlet_convergence.png'));
hold off;
