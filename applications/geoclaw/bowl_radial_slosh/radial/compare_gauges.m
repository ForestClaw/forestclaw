close all;

mq = 1;    % 1=h; 2=hu; 3=hv; 4=eta
for n = 1:9,
    plot_gauge(n,mq);
    hold on;
    h = plot_gauge(n+100);
    set(h,'color','r');
    shg;
    input('Hit enter to see next gauge : ');
    hold off;
end