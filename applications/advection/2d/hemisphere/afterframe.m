setviews;

% showgridlines(1:4)
showpatchborders();

yrbcolormap;

view(3)
daspect([1 1 1]);
axis off;

axis([-1.1 1.1 -1.1 1.1 0 1.1]);

NoQuery = 0;
prt = false;
if (prt)
  filename = sprintf('hemisphere%03d.png',Frame);
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;


clear afterframe;
clear mapc2m;
clear parallelpartitions
