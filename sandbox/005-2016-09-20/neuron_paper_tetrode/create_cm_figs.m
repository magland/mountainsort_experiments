close all;

font_size=14;
fs2=20;

conf_count_plot(csvread('manual1-vs-manual2.csv'),font_size,'');
xlabel('Manual 2','FontSize',fs2);
ylabel('Manual 1','FontSize',fs2);
saveas(gcf,'manual1-vs-manual2.png');

conf_count_plot(csvread('manual1-vs-manual3.csv'),font_size,'');
xlabel('Manual 3','FontSize',fs2);
ylabel('Manual 1','FontSize',fs2);
saveas(gcf,'manual1-vs-manual3.png');

conf_count_plot(csvread('manual2-vs-manual3.csv'),font_size,'');
xlabel('Manual 3','FontSize',fs2);
ylabel('Manual 2','FontSize',fs2);
saveas(gcf,'manual2-vs-manual3.png');

conf_count_plot(csvread('manual1-vs-ms.csv'),font_size,'');
xlabel('MountainSort','FontSize',fs2);
ylabel('Manual 1','FontSize',fs2);
saveas(gcf,'manual1-vs-ms.png');

conf_count_plot(csvread('manual2-vs-ms.csv'),font_size,'');
xlabel('MountainSort','FontSize',fs2);
ylabel('Manual 2','FontSize',fs2);
saveas(gcf,'manual2-vs-ms.png');

conf_count_plot(csvread('manual3-vs-ms.csv'),font_size,'');
xlabel('MountainSort','FontSize',fs2);
ylabel('Manual 3','FontSize',fs2);
saveas(gcf,'manual3-vs-ms.png');