s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders(1:10);
setpatchborderprops('linewidth',1);
caxis([-1,1])

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
