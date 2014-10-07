s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders(1:10);
caxis([0,1])

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end;

shg
clear afterframe
clear mapc2m
