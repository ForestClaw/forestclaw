global iscart;

if (iscart)
    s = 1e-2;
else
    s = alpha+1e-2;
end
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders(1:10);
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
daspect([1,1,1]);

if (iscart)
    view(2);
else
    view(3);
end

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
clear mapc2m_torus;
