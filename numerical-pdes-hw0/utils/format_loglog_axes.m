function format_loglog_axes()
%FORMAT_LOGLOG_AXES Apply consistent formatting to log-log plots.
%
% Sets LaTeX tick labels, increases font size, enables grid lines, and
% ensures consisten styling across convergence plots.
grid on;
ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.TickLabelInterpreter = 'latex';
end
