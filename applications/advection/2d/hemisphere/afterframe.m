setviews;
s = 0.1;
axis([-1-s 1+s -1-s 1+s 0 1+s]);

showpatchborders();
setpatchborderprops('linewidth',1);

yrbcolormap;
caxis([0 1]);

view(3)
daspect([1 1 1]);
axis off;


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
