function h = report(param,out)
%REPORT(PARAM,OUT) Show a simple report.

LineWidth = 2;
h = figure;
subplot(2,1,1);

if isfield(out,'err'),
    semilogy(param.restart_length*(1:length(out.err)),out.err,'k-o','LineWidth',LineWidth);
    hold on
end;
if isfield(out,{'bound_1','bound_2'}),
    semilogy(param.restart_length*(1:length(out.bound_1)),out.bound_1,'g:o','LineWidth',LineWidth);
    semilogy(param.restart_length*(1:length(out.bound_2)),out.bound_2,'g:o','LineWidth',LineWidth);
    hold on
end;
if isfield(out,'update'),
    semilogy(param.restart_length*(1:length(out.update)),out.update,'r:o','LineWidth',LineWidth);
    hold on
end;
grid on
xlabel('mvp');
ylabel('absolute error');
title('error (black), estimators (green-dot), norm of updates (red-dot)');

%%%%%%

subplot(2,1,2);
if isfield(out,'time'),
    plot(param.restart_length*(1:length(out.time)),out.time,'c-o','LineWidth',LineWidth);
    hold on
    plot(param.restart_length*(1:length(out.time)),cumsum(out.time),'b-o','LineWidth',LineWidth);
    xlabel('mvp');
    ylabel('seconds');
    title('execution time per cycle (cyan) and cumulative (blue)');
    x = get(gca,'XLim'); x = x(1);
    y = get(gca,'YLim'); y = y(end);
    text(x,0.9*y,[' Stop condition: ' out.stop_condition],'Interpreter','none');
end;