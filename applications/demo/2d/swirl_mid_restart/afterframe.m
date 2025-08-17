s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders;
setpatchborderprops('linewidth',1)

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


view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg
