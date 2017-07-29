function plot_velocity(N)

close all;

% Quiver plot
figure(1);

x = linspace(0,1,N+1);
y = x;
[xm,ym] = meshgrid(x,y);

[um,vm] = vel(xm,ym);

s = 1;
quiver(xm,ym,um,vm,s);
daspect([1 1 1]);
axis([0 1 0 1]);
title('Vecocity','fontsize',18);

% Contour plot of magnitude
figure(2);
clf

x = linspace(0,1,10*N+1);
y = x;

[xm,ym] = meshgrid(x,y);

[um,vm] = vel(xm,ym);

mm = sqrt(um.^2 + vm.^2);
contour(xm,ym,mm,11);
hold on;
contour(xm,ym,mm,[0 0],'r--');
daspect([1 1 1]);
axis([0 1 0 1]);
colorbar;

title('Speed','fontsize',18);


% Divergence plot
figure(3)
% plot divergence
xe = linspace(0,1,10*N + 1);
xc = (xe(1:end-1) + xe(2:end))/2;

ye = linspace(0,1,10*N + 1);
yc = (ye(1:end-1) + ye(2:end))/2;

[xem,ycm] = meshgrid(xe,yc);
[um,~] = vel(xem,ycm);

[xcm,yem] = meshgrid(xc,ye);
[~,vm] = vel(xcm,yem);

h = xe(2) - xe(1);

ux = (um(:,2:end) - um(:,1:end-1))/h;
vy = (vm(2:end,:) - vm(1:end-1,:))/h;

divu = ux + vy;

[xcm,ycm]= meshgrid(xc,yc);

contour(xcm,ycm,divu,10);
daspect([1 1 1]);
colorbar;
title('Divergence','fontsize',18);


end

function [u,v] = vel(x,y)


ucc = @(x,y) 2*((sin(pi*x)).^2 .* sin(pi*y) .* cos(pi*y));
vcc = @(x,y) -2*sin(pi*y).^2 .* sin(pi*x) .* cos(pi*x);
vcc_add = @(x,y) 0.2*sin(4*pi*y).*cos(4*pi*x);

ucc = @(x,y) cos(pi*x).^2 - 0.5;
vcc = @(x,y) sin(pi*y).^2 - 0.5;

u = ucc(x,y);
v = vcc(x,y);


end
