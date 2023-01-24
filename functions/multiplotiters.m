function multiplotiters(PLOTS1,PLOTS2)
figure;
hold on;
plot(PLOTS1.ITERS, PLOTS1.FVALS,'-rs', 'LineWidth', 2,...
    'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'g',...
    'MarkerSize', 5);
plot(PLOTS2.ITERS, PLOTS2.FVALS,':bo', 'LineWidth', 2,...
    'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'r',...
    'MarkerSize', 5);
xlabel('Iteration'),ylabel('Objectif');
legend('DCA','Local DCA');
hold off;
end