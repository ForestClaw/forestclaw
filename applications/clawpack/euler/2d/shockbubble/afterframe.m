if (PlotType == 3)
    caxis([0 200]);
elseif mq == 5
    % 0=y; 1=b; g=2; r=3;
    b = [0 0 1];
    g = [0 0.8 0];
    r = [1 0 0];
    y = [1 1 0];
    s = linspace(0,1,11)';
    yb = kron(y,1-s) + kron(b,s);
    bg = kron(b,1-s) + kron(g,s);
    gr = kron(g,1-s) + kron(r,s);
    ry = kron(r,1-s) + kron(y,s);
    cm = [y; yb; b; bg; g; gr; r];
    cm = [r; ry; y; yb; b; bg; g];
    colormap((cm));
    
    
    caxis([0.1 3.5])
    % yrbcolormap;
else
    colormap(jet)
    % caxis([0.1 2.81]);
    caxis([0.1 3.5])
end

% Axes and title
axis([0 2 0 0.5])
set(gca,'fontsize',16);
tstr = sprintf('ForestClaw : t = %12.4f',t);
title(tstr,'fontsize',16);

% Show patch borders
showpatchborders
setpatchborderprops('linewidth',1);

% This is used for creating vectorized PDFs
prt_tikz = false;
if (prt_tikz)
    axis off
    hidepatchborders;
    delete(get(gca,'title'));
    figsize = [4,1];  % Should match size set in options
    set(gcf,'papersize',figsize);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 figsize]);
    fname = sprintf('plot_%04d.png',Frame);
    
    % Match print resolution to computational resolution
    maxlevel = 4;   % Adjust to determine resolution
    res = sprintf('-r%d',mx*2^maxlevel);
    print('-dpng',res,fname);  
end
daspect([1,1,1]);
shg

clear afterframe
clear mapc2m
clear parallelpartitions
