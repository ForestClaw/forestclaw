function plot_cyl_in_sphere()

clf;

N = 200;

x = linspace(0,1,N);
y = x;
[xm,ym] = meshgrid(x,y);

% Sphere
R = sqrt(pi^2 + 1);
arc = pi/2 - acos(pi/R);
phi0 = -0.5*arc/(pi/2);
phi1 = -phi0;

[xp,yp,zp] = mapc2m_latlong2(xm,ym,R,phi0,phi1);
p = patch(surf2patch(xp,yp,zp));
set(p,'edgecolor','none','facecolor','b','facealpha',0.25);

hold on;

% Cylinder
[xp,yp,zp] = mapc2m(xm,ym);
p = patch(surf2patch(xp,yp,zp-pi));
set(p,'edgecolor','none','facecolor','y');

set(gca,'clipping','off');
daspect([1,1,1]);

view(3);

camlight;


shg;



end