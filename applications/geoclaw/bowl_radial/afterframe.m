s = 1e-2;
axis([-100 100 -100 100])
daspect([1 1 1]);
axis on;

showpatchborders(5:10);
set(gca,'zlim',[-8,0]);   % Need so that all patchborders show up

caxis([-0.9,0.9])

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'bowl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
