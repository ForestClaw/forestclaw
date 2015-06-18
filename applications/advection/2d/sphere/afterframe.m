setviews;
daspect([1 1 1]);
axis off;


yrbcolormap;
caxis([0 1]);
showgridlines(1:2);
showpatchborders(1:7);
if (ShowUnderOverShoots == 1)
    under_label = sprintf('%6.2e',qmin);
    over_label = sprintf('1 + %6.2e',qmax-1);
    colorbar_underover(under_label,over_label);
end

zoom(1.1^5);
hidegridlines;

NoQuery = 0;
prt = false;
if (prt)
    filename = sprintf('sphere%4.4d.png',Frame);
      disp(filename);
      print('-dpng',filename);
end;

shg;


clear afterframe;
clear mapc2m;
clear mapc2m_cubed_sphere;
clear parallelpartitions;
