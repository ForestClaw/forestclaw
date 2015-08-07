s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

colormap(jet);
cv = linspace(qmin,qmax,21);
cv([1 end]) = [];
% drawcontourlines(cv);
% setcontourlineprops('linewidth',2);
setpatchborderprops('linewidth',1);
showpatchborders;

% caxis([-1e-8,1e-8]);
% colormap(yrbcolormap);

fprintf('%10s : %12.4e\n','qmin',qmin);
fprintf('%10s : %12.4e\n','qmax',qmax);

view(2);

NoQuery = 0;
prt = false;
if (prt)
  str = sprintf('quadrant%02d.png',Frame);
  fprintf('Printing file %s\n',str);
  print('-dpng',str);
end;

shg

clear afterframe;
clear parallelpartitions;
