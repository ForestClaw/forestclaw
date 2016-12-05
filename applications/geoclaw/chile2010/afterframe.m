% chile parameters

axis([-120 -60 -60 0])
daspect([1 1 1]);
axis on;

m = min(abs([qmin,qmax]));
caxis([-0.2 0.2]);

showpatchborders;

hold on;
add_gauges();

hold on;
plot_regions(t);
hold off;

set(gca,'fontsize',16);

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'chile0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
