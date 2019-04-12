s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

showpatchborders;

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
clear mapc2m_fivepatch;
