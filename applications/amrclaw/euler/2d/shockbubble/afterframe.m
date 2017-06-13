if (PlotType ~= 3)
    colormap(jet)
end

%%
axis([0 2 0 0.5])


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
    
    
    caxis([0 3])
    % yrbcolormap;
else
    colormap(jet)
    caxis([0.1 2.81]);
end

showpatchborders
setpatchborderprops('linewidth',1);
% showgridlines(1:5);
% hidepatchborders;
set(gca,'fontsize',16);

tstr = sprintf('ForestClaw : t = %12.4f',t);
title(tstr,'fontsize',16);

showpatchborders;
% colormap(white)


shg;

clear afterframe
clear mapc2m
clear parallelpartitions
