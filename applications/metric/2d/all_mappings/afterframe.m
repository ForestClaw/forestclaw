axis image;
daspect([1 1 1]);

colormap(jet);

hidepatchborders;
setviews;
axis off;

showpatchborders;
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


fprintf('%20s %16.8e\n','qmin',qmin);
fprintf('%20s %16.8e\n','qmax',qmax);

NoQuery = 0;
prt = false;
if (prt)
  filename = 'sphere000.png';
  str = num2str(Frame);
  len = length(str);
  filename(8-len+1:8) = str;
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;

shg;

clear afterframe;
clear mapc2m;
clear underover;
