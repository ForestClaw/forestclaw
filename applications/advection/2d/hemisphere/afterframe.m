s = 1e-2;
% axis([-1-s 1+s -1-s 1+s -1-s 1+s])
axis image;
daspect([1 1 1]);
axis off;

showgridlines(1:4)
setpatchborderprops(1:6,'linewidth',2);
showpatchborders();
setviews;

yrbcolormap;
if (ShowUnderOverShoots == 1)
  under_label = sprintf('%6.2e',qmin);
  over_label = sprintf('1 + %6.2e',qmax-1);
  underover_colorbar(under_label,over_label);
  fprintf('%12s %15s\n','qmin',under_label);
  fprintf('%12s %15s\n\n','qmax',over_label);
else
  fprintf('%12s %24.16f\n','qmin',qmin);
  fprintf('%12s %24.16f\n\n','qmax',qmax);
end;



NoQuery = 0;
prt = false;
if (prt)
  filename = 'swirl000.png';
  str = num2str(Frame);
  len = length(str);
  filename(8-len+1:8) = str;
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
