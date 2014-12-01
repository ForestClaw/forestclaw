s = 0;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);

if (PlotParallelPartitions == 0)
    yrbcolormap;
end
showpatchborders(1:9);
setpatchborderprops('linewidth',2);
% hidepatchborders;
% delete(get(gca,'title'));
caxis([0 1]);

% colormap(white);

caxis([0,1])
qlo = 0;
qhi = 1;
under_label = sprintf('0 - %7.1e',qlo-qmin);
over_label = sprintf('1 + %7.1e',qmax-qhi);
fprintf('%6s %12s\n','qmin',under_label);
fprintf('%6s %12s\n\n','qmax',over_label);

if (ShowUnderOverShoots)
    qlo = 0;
    qhi = 1;
    colorbar_underover(under_label,over_label);
end


shg;

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'filament0000','png');    
    print('-dpng','-r1600',filename);
end;

clear afterframe
clear mapc2m
clear mapc2m_squareddisk
clear mapc2m_pillowdisk

