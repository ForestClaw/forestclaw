s = 1e-2;
axis([-100 100 -100 100])
daspect([1 1 1]);
axis on;

showpatchborders(1:6);
hidepatchborders(7)
cv = linspace(qmin,qmax,20);
% drawcontourlines(cv);
set(gca,'zlim',[-8,0]);   % Need so that all patchborders show up

rybcolormap;
tol = 1e-1;
caxis([-tol,tol])
% caxis([qmin,qmax])
daspect([20,20,1]);

hold on;
zmin = 80;
ep = 0.01;
w = sqrt(zmin/ep);
th = linspace(0,2*pi,500);
plot(w*cos(th),w*sin(th),'k','linewidth',2);
hold off;

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
