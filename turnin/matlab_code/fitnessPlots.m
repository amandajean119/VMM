filename1 = '~/Git/VMM/src/fitness0mut0crossnoElitism';
filename2 = '~/Git/VMM/src/fitness08mut0crossnoElitism';
filename3 = '~/Git/VMM/src/fitness0mut08crossnoElitism';
filename4 = '~/Git/VMM/src/fitness08mut08crossnoElitism';
filename5 = '~/Git/VMM/src/fitness0mut0crossElitism';
filename6 = '~/Git/VMM/src/fitness08mut0crossElitism';
filename7 = '~/Git/VMM/src/fitness0mut08crossElitism';
filename8 = '~/Git/VMM/src/fitness08mut08crossElitism';
csv1 = csvread(strcat(filename1, '.txt'));
csv2 = csvread(strcat(filename2, '.txt'));
csv3 = csvread(strcat(filename3, '.txt'));
csv4 = csvread(strcat(filename4, '.txt'));
csv5 = csvread(strcat(filename5, '.txt'));
csv6 = csvread(strcat(filename6, '.txt'));
csv7 = csvread(strcat(filename7, '.txt'));
csv8 = csvread(strcat(filename8, '.txt'));


hold on;
plot(1:length(csv1(:,3)), csv1(:,3), 'b', 'LineWidth', 1.2);
plot(1:length(csv1(:,4)), csv1(:,4), '--b', 'LineWidth', 1.2);
plot(1:length(csv5(:,3)), csv5(:,3), 'g', 'LineWidth', 1.2);
plot(1:length(csv5(:,4)), csv5(:,4), '--g', 'LineWidth', 1.2);
set(gca, 'fontsize', 20)
title('0.1 Crossover rate, 0.1 mutation rate ', 'fontsize', 25);
axis([0, 50, 0.0, 1.0]);
xlabel('Generations','FontSize', 25);
ylabel('Fitness','FontSize', 25);
h_legend = legend('No elitism Best Fitness', 'No elitism Avg Fitness', 'With elitism Best Fitness', 'With elitism Avg Fitness','Location','SouthEast');
set(h_legend,'FontSize',15);
print('-dpng', '-r200', filename1);
clf;

hold on;
plot(1:length(csv2(:,3)), csv2(:,3), 'b', 'LineWidth', 1.2);
plot(1:length(csv2(:,4)), csv2(:,4), '--b', 'LineWidth', 1.2);
plot(1:length(csv6(:,3)), csv6(:,3), 'g', 'LineWidth', 1.2);
plot(1:length(csv6(:,4)), csv6(:,4), '--g', 'LineWidth', 1.2);
set(gca, 'fontsize', 20)
title('0.1 Crossover rate, 0.8 mutation rate  ', 'fontsize', 25);
axis([0, 50, 0.0, 1.0]);
xlabel('Generations','FontSize', 25);
ylabel('Fitness','FontSize', 25);
h_legend = legend('No elitism Best Fitness', 'No elitism Avg Fitness', 'With elitism Best Fitness', 'With elitism Avg Fitness','Location','SouthEast');
set(h_legend,'FontSize',15);
print('-dpng', '-r200', filename2);
clf;

hold on;
plot(1:length(csv3(:,3)), csv3(:,3), 'b', 'LineWidth', 1.2);
plot(1:length(csv3(:,4)), csv3(:,4), '--b', 'LineWidth', 1.2);
plot(1:length(csv7(:,3)), csv7(:,3), 'g', 'LineWidth', 1.2);
plot(1:length(csv7(:,4)), csv7(:,4), '--g', 'LineWidth', 1.2);
set(gca, 'fontsize', 20)
title('0.8 Crossover rate, 0.1 mutation rate  ', 'fontsize', 25);
axis([0, 50, 0.0, 1.0]);
xlabel('Generations','FontSize', 25);
ylabel('Fitness','FontSize', 25);
h_legend = legend('No elitism Best Fitness', 'No elitism Avg Fitness', 'With elitism Best Fitness', 'With elitism Avg Fitness','Location','SouthEast');
set(h_legend,'FontSize',15);
print('-dpng', '-r200', filename3);
clf;

hold on;
plot(1:length(csv4(:,3)), csv4(:,3), 'b', 'LineWidth', 1.2);
plot(1:length(csv4(:,4)), csv4(:,4), '--b', 'LineWidth', 1.2);
plot(1:length(csv8(:,3)), csv8(:,3), 'g', 'LineWidth', 1.2);
plot(1:length(csv8(:,4)), csv8(:,4), '--g', 'LineWidth', 1.2);
set(gca, 'fontsize', 20)
title('0.8 Crossover rate, 0.8 mutation rate  ', 'fontsize', 25);
axis([0, 50, 0.0, 1.0]);
xlabel('Generations','FontSize', 25);
ylabel('Fitness','FontSize', 25);
h_legend = legend('No elitism Best Fitness', 'No elitism Avg Fitness', 'With elitism Best Fitness', 'With elitism Avg Fitness','Location','SouthEast');
set(h_legend,'FontSize',15);
print('-dpng', '-r200', filename4);
clf;