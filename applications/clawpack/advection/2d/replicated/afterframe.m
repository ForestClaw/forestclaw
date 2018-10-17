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

plot_tikz = true;
if (plot_tikz)
    caxis([0 7]);
    axis([0, 1, 0,1]);
    figsize = [8,8];  % Should match tikz figsize.
    maxlevel = 7;
    dpi = mx*2^maxlevel/figsize(1);
    prefix = 'plot';
    caxis([0,1]);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end

shg

clear afterframe;
clear mapc2m;
