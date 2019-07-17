
% Axes and title
axis([0 2 0 0.5])
set(gca,'fontsize',16);
tstr = sprintf('ForestClaw : t = %12.4f',t);
title(tstr,'fontsize',16);
daspect([1,1,1]);

% Color map and axis
colormap(jet)
caxis([0, 4])

% Show patch borders
showpatchborders
setpatchborderprops('linewidth',1);

% This is used for creating vectorized PDFs
prt_tikz = false;
if (prt_tikz)
    hidepatchborders;
    mi = 4;
    mj = 1;
    figsize = [4,1];  % Should match tikz figsize.
    maxlevel = 6;
    dpi = mi*mx*2^maxlevel/figsize(1);
    prefix = 'plot';
    caxis([0,4]);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end
shg
