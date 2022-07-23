% Solve the filament problem using Lagrangian method. 

function [xout,yout] = filament_soln(N,T)

global N_filament dir_filament;
N_filament = N;

reverse_flow = false;

th = linspace(0,2*pi,N)';

R_init = 0.25;
xc_init = 0.5;
yc_init = 1;

% Initial disk of radius
x = R_init*cos(th) + xc_init;
y = R_init*sin(th) + yc_init;

dir = 1;
dir_filament = dir;

p0 = [x;y];
% opt.AbsTol = 1e-12;
opt.RelTol = 1e-12;
[tout,pout] = ode45(@f_rhs,[0,T],p0,opt);

xout = pout(end,1:N);
yout = pout(end,N+1:end);
if (nargout > 0)
    return;
end


% plot forward flow
figure(2);
h = fill(xout,yout,'r');
set(h,'edgecolor','none');
hold on;
plot(xout,yout,'k-','linewidth',1);
axis([0 2 0 2]);
daspect([1 1 1]);

% Area of the initial disk of radius 0.25, centered at (0.5,1).  
area_init = pi*R_init^2;
area_init2 = polyarea(x,y);

area_end = polyarea(xout,yout);
fprintf('Forward flow\n');
fprintf('%20s %24.16f\n','Initial area',area_init);
fprintf('%20s %24.16f\n','Initial area',area_init2);
fprintf('%20s %24.16f\n','Final area',area_end);
fprintf('%20s %24.4e\n','Error',abs(area_init2 - area_end));
set(gca,'fontsize',16);

if (reverse_flow)
    % Reverse the flow
    dir = -1;
    figure(3);
    p0 = pout(end,:)';
    area_init = pi*R_init^2;
    [tout,pout] = ode45(@f_rhs,[0,T],p0,opt);

    xout = pout(end,1:N);
    yout = pout(end,N+1:end);
    h = fill(xout,yout,'r');
    set(h,'edgecolor','none');
    hold on;
    plot(xout,yout,'k-','linewidth',1);
    axis([0 2 0 2]);
    daspect([1 1 1]);
    hold on;
    plot(x,y,'k-','linewidth',2);
end
    
area_init = polyarea(xout,yout);
fprintf('\n');
fprintf('Backward flow\n');
fprintf('%20s %24.16f\n','Initial area',area_init);
fprintf('%20s %24.16f\n','Initial area',area_init2);
fprintf('%20s %24.16f\n','Final area',area_end);
fprintf('%20s %24.4e\n','Error',abs(area_init2 - area_end));
set(gca,'fontsize',16);

end

function f = f_rhs(t,p)
    
global N_filament dir_filament;
N = N_filament;
dir = dir_filament;

x = p(1:N);
y = p(N+1:end);
% (4.d0/3.d0)*r**3
r = sqrt((x-1).^2 + (y-1).^2);
rx = (x-1)./r;
ry = (y-1)./r;

% u = dPsi/dy  = d(4/3*r^3)/dr*drdy
% v = -dPsi/dx = -d(4/3*r^3)/dr*drdx
u = 4*r.^2.*ry;
v = -4*r.^2.*rx;

f = dir*[u; v];

end



