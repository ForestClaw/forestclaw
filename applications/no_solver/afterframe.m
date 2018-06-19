s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders(1:10);
setpatchborderprops('linewidth',1);
caxis([0,1]) 

plot_tikz = true;
if (plot_tikz)
    axis([-1, 1, -1,1]);
    figsize = [16,16];  % Should match tikz figsize.
    maxlevel = 6;
    dpi = mx*2^maxlevel/figsize(1);
    fprintf('dpi = %d\n',dpi);
    prefix = 'plot';
    caxis([0,4]);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'mesh0000','png');
  print('-dpng',filename);
end

shg
clear afterframe
clear mapc2m
