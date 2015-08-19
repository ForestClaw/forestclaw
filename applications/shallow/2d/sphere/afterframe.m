setviews;

global map isflat;

alpha = 0.4;
s = 1e-2;    
alim = [-1 1];
alim = alim + [-s s];
axis([alim alim]);
daspect([1 1 1]);

if (PlotParallelPartitions == 0)
    caxis([0.2 0.21]);
end

yrbcolormap;
showpatchborders(1:10);

fprintf('%10s : %12.4e\n','qmin',qmin);
fprintf('%10s : %12.4e\n','qmax',qmax);

colorbar;

daspect([1,1,1]);

if Frame >= 18
    view(vbot);
else
    view(3);
    zoom(1.5);
end

axis off
% colorbar

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'torus0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
clear parallelpartitions;
