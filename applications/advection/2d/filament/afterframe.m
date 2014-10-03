s = 0;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);

if (PlotParallelPartitions == 0)
    yrbcolormap;
end
showpatchborders(1:8);
% hidepatchborders;
delete(get(gca,'title'));
caxis([0 1]);

% colormap(white);


shg;

NoQuery = 0;
prt = false;
if (prt)
    MaxFrames = 65;
    filename = framename(Frame,'filament0000','png');
    print('-dpng',filename);
end;

clear afterframe
clear mapc2m
clear mapc2m_squareddisk
clear mapc2m_pillowdisk

