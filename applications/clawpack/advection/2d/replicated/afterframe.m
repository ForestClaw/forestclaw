yrbcolormap;
setviews;

global map isflat;

brick = load('brick.dat');
axis([0 brick(1,1) 0 brick(1,2)])
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1)

view(2);
delete(get(gca,'title'));
set(gca,'fontsize',16,'box','on');
axis square

% This is used for creating vectorized PDFs
prt_tikz = true;
if (prt_tikz)
    axis off
    hidepatchborders;
    delete(get(gca,'title'));
    figsize = [1,1];  % Should match size set in options
    set(gcf,'papersize',figsize);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 figsize]);
    fname = sprintf('plot_%04d.png',Frame);
    
    % Match print resolution to computational resolution
    maxlevel = 7;   % Adjust to determine resolution
    res = sprintf('-r%d',mx*2^maxlevel);
    print('-dpng',res,fname);  
end


shg

clear afterframe;
clear mapc2m;
