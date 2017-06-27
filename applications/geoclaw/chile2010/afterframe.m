setviews;

% Patches
showpatchborders;
setpatchborderprops('linewidth',1);

% Color maps, etc
if (PlotParallelPartitions)
    if (Frame == 0)
        cm = ppcolors(4);
    end
    colormap(cm);
    ppcolors_colorbar(4);
else
    colormap(parula);
    colorbar;    
end
caxis([-0.2 0.2]);

% Gauges, regions
hold on;
add_gauges();
add_regions(t);
hold off;

fprintf('\n');
fprintf('qmin : %12.4e\n',qmin);
fprintf('qmax : %12.4e\n',qmax);
fprintf('\n');


% Axes
axis([-120 -60 -60 0])
daspect([1 1 1]);
title(sprintf('Surface at time %6.2f hours',t/(60^2)),'fontsize',18);
set(gca,'fontsize',16);
axis on;


NoQuery = 0;
prt = false;
if (prt)
    filename = sprintf('chile%04d.png',Frame);
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
