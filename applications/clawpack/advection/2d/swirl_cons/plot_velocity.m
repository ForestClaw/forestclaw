function plot_velocity(ex,M)

global example

close all;

N = 200;  % For hi-res plots

example = ex;

% Quiver plot
figure(1);

x = linspace(0,1,M+1);
y = x;
[xm,ym] = meshgrid(x,y);

[um,vm] = vel(xm,ym);

s = 1.5;
quiver(xm,ym,um,vm,s);
daspect([1 1 1]);
axis([0 1 0 1]);
title(sprintf('Example %d : Velocity',ex),'fontsize',18);
set(gca,'fontsize',16);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);

% Contour plot of magnitude
figure(2);
clf

x = linspace(0,1,N+1);
y = x;

[xm,ym] = meshgrid(x,y);

[um,vm] = vel(xm,ym);

mm = sqrt(um.^2 + vm.^2);
contour(xm,ym,mm,21);

fprintf('%-10s %12.4e\n','min(mm)',min(mm(:)));
fprintf('%-10s %12.4e\n','max(mm)',max(mm(:)));


caxis([0 1]);
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
dmax = round(max(divu));
dmin = round(min(divu));
dlim = max(abs([dmax,dmin]));

[xcm,ycm]= meshgrid(xc,yc);

contour(xcm,ycm,divu,21);
daspect([1 1 1]);
caxis([-dlim dlim]);
colorbar;
title('Divergence','fontsize',18);


end

function [u,v] = vel(x,y)

global example

if (example == 1)
    a = 0.5;
else
    a = -0.5;
end

s = sqrt(2)/2;
ucc = @(x,y) s*(cos(pi*x).^2 + a);
vcc = @(x,y) s*(sin(pi*y).^2 + a);

u = ucc(x,y);
v = vcc(x,y);


end
