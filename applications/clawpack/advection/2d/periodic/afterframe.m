s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders(1:10);
setpatchborderprops('linewidth',1);
caxis([0,1])
view(2);

NoQuery = 0;
prt = false;
if (prt)
    hidepatchborders;
    dpi = 2^7;
    figsize = [4,4];
    prefix = 'plot';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
end

shg

clear afterframe
