% chile parameters

axis([-120 -60 -60 0])
daspect([1 1 1]);
axis on;

m = min(abs([qmin,qmax]));
caxis([-0.2 0.2]);

showpatchborders;

hold on;
add_gauges();
add_regions(t);
hold off;

fprintf('\n');
fprintf('qmin : %12.4e\n',qmin);
fprintf('qmax : %12.4e\n',qmax);
fprintf('\n');

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
