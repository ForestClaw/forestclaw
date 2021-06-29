function plot_cyl_sphere(xbar)

close all;

setviews;

R_cyl = 1;
H_cyl = 2*pi;


% Mapping parameters
if (nargin == 0)
    xbar = R_cyl;
end

% Create mesh
mx = 70;
xe = linspace(0,1,mx+1);

my = 20;
ye = linspace(0,1,my+1);

[xem,yem] = meshgrid(xe,ye);

% ----------------------------------
% Mapping
% ----------------------------------

% Check that we get a true cylinder if xbar == R_cyl : 

[x0,y0,z0] = mapc2m_cyl_sphere(0,0,'cylinder',R_cyl,H_cyl);
[x1,y1,z2] = mapc2m_cyl_sphere(0,0.5,'cylinder',R_cyl,H_cyl);
[x2,y2,z2] = mapc2m_cyl_sphere(0,1,'cylinder',R_cyl,H_cyl);

r2 = sqrt(x2^2 + y2^2);
fprintf('%20s %24.16f\n','r (top)',r2);

r1 = sqrt(x1^2 + y1^2);
fprintf('%20s %24.16f\n','r (middle)',r1);

r0 = sqrt(x0^2 + y0^2);
fprintf('%20s %24.16f\n','r (bottom edge)',r0);


if (abs(r0-r1)/r1 > eps(1))
    error('Not a true cylinder');
end


% % Plot cylinder
[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,'cylinder',R_cyl,H_cyl);

p = patch(surf2patch(xep,yep,zep));
set(p,'edgecolor','none','facecolor','r','facealpha',0.1);
hold on;

% plot top and bottom of cylinder
th = linspace(0,2*pi,500);
plot3(R_cyl*cos(th),R_cyl*sin(th),H_cyl/2 + 0*th,'r','linewidth',2);
hold on;
plot3(R_cyl*cos(th),R_cyl*sin(th),-H_cyl/2 + 0*th,'r','linewidth',2);

hold on;

% Plot sphere
[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,'sphere',R_cyl,H_cyl);

p = patch(surf2patch(xep,yep,zep));
set(p,'edgecolor','none','facecolor','b','facealpha',0.15);


% User shape
[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,xbar,R_cyl,H_cyl);

p = patch(surf2patch(xep,yep,zep));
set(p,'edgecolor','k','facecolor','none','facealpha',0.15);

% ----------------------------------
% Misc axis items
% ----------------------------------
set(gca,'clipping','off');
daspect([1,1,1]);
% camlight;

axis off;
view(3);

set(gcf,'position',[1649,386, 560, 420]);
shg



end