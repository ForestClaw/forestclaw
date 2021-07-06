function plot_cyl_in_sphere()

R_cyl = 1;
H_cyl = 2*pi*R_cyl;

% Sphere
R_latlong = sqrt((H_cyl/2)^2 + R_cyl^2);
phi_circle = asin(H_cyl/(2*R_latlong));

% Map 0 --> -phi_circle/(pi);   
%     1 -->  phi_circle/(pi);  
%    l(y) = 2*(phi_circle/pi)*y - phi_circle/(pi) = phi_circle/(pi)*(2*y - 1)
phi0 = -phi_circle/pi;
phi1 = -phi0;

% ---------------------------- Plotting -----------------------------------
figure(1);
clf;

set_blocknumber(0);

N = 200;

x = linspace(0,1,N);
y = x;
[xm,ym] = meshgrid(x,y);

set_blocknumber(0);
R_cyl = 1;
H_cyl = 2*pi*R_cyl;

% Sphere
R_latlong = sqrt((H_cyl/2)^2 + R_cyl^2);
arc = asin(H_cyl/(2*R_latlong));
phi0 = -0.5*arc/(pi/2);
phi1 = -phi0;

[xp,yp,zp] = mapc2m_latlong2(xm,ym,R_latlong,phi0,phi1);
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