s = 1e-2;
axis([-100 100 -100 100])
daspect([1 1 1]);
axis on;

showpatchborders(3:10);
hidepatchborders(7)
cv = linspace(qmin,qmax,20);
% drawcontourlines(cv);
set(gca,'zlim',[-8,0]);   % Need so that all patchborders show up

caxis([-0.9,0.9])

view(2);
colorbar
NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'bowl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
