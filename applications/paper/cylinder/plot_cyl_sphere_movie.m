function plot_cyl_sphere_movie()

close all;

R_cyl = 2;
H_cyl = 2*pi;

mx = 70;
xe = linspace(0,1,mx+1);

my = 20;
ye = linspace(0,1,my+1);

[xem,yem] = meshgrid(xe,ye);

[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,'cylinder',R_cyl,H_cyl);

% Cylinder
q(1) = patch(surf2patch(xep,yep,zep));
set(q(1),'edgecolor','none','facecolor','r','facealpha',0.2);

hold on;

% Sphere
[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,'sphere',R_cyl,H_cyl);

q(2) = patch(surf2patch(xep,yep,zep));
set(q(2),'edgecolor','none','facecolor','b','facealpha',0.2);


axis off;

daspect([1,1,1]);

set(gcf,'position',[1649,386, 560, 420]);
s = 8;

set(gca,'clipping','off');

axis([-s,s,-s,s]);

setviews;
view([-37.322436719317636,25.641518552290719]);




% create sequence for xbar.  Negative values aren't really doing the right
% thing yet.  
xbar_vec = linspace(-10,10,61);

[xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,xbar_vec(1),R_cyl,H_cyl);

FVC = surf2patch(xep,yep,zep);
p = patch(FVC);
set(p,'edgecolor','k','facecolor','none','facealpha',0.1);
hold on;

shg;


for i = 1:length(xbar_vec)
    xbar = xbar_vec(i);
    [xep,yep,zep] = mapc2m_cyl_sphere(xem,yem,xbar,R_cyl,H_cyl);

    FVC = surf2patch(xep,yep,zep);
    set(p,'Faces',FVC.faces,'vertices',FVC.vertices);
    
    axis([-s,s,-s,s]);
    view([-37.322436719317636,25.641518552290719]);
    pause(0.1);
end
   
end