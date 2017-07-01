% Patch borders
showpatchborders;
setpatchborderprops('linewidth',1);

% Contour lines
cv = linspace(qmin,qmax,11);
%drawcontourlines(cv);
%setcontourlineprops('color','r');

% Color maps
colormap(parula)
caxis([-1 2])
colorbar;

% Axes
axis([-100 100 -100 100])
set(gca,'zlim',[-8,0]);   % Need so that all patchborders show up
daspect([1 1 1]);
set(gca,'fontsize',16);
axis on;

NoQuery = 0;
prt = false;
if (prt)
  filename = sprintf('flat%04d.png',Frame);
  print('-dpng',filename);
end

shg

clear afterframe;
