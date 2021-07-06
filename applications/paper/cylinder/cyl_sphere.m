function cyl_sphere(xbar)

close all;

R_cyl = 1;
H_cyl = 2*pi;

beta = H_cyl/R_cyl;
if (nargin == 0)
    xbar = R_cyl*sqrt(1 + beta^2/4);
end

if (xbar == R_cyl)
    c = beta^2/4;
    e = sqrt(eps(c)/10);
    xbar = R_cyl*(1 + e);
end

alpha = xbar/R_cyl;

% (-L,0) is the center of the sphere.
L        = -R_cyl/2*(1 + alpha + beta^2/(4*(1-alpha)));
R_sphere = -R_cyl/2*(1 - alpha + beta^2/(4*(1-alpha)));

% Get range for phi. 
% TODO : Fix so that angle flips at L=0.  
if (-L <= R_cyl)
    phi0 = asin(-H_cyl/(2*R_sphere));
    phi1 = -phi0;
else
    phi1 = acos((L+R_cyl)/R_sphere);
    phi0 = -phi1;
end


fprintf('\n');
fprintf('%15s  %24.16f\n','min xbar',1);
fprintf('%15s  %24.16f\n','xbar',xbar);
fprintf('%15s  %24.16f\n','max xbar',sqrt(1 + beta^2/4));
fprintf('\n');
fprintf('%15s  %24.16f\n','phi0',phi0);
fprintf('%15s  %24.16f\n','phi1',phi1);
fprintf('\n');


fprintf('%15s  %24.17g\n','R_sphere',R_sphere);
fprintf('%15s  %24.17g\n','-L (xcenter)',-L);

phi = linspace(phi0,phi1,501);

xc = -L;
yc = 0;

x = R_sphere*cos(phi) + xc;   % This is the radius
y = R_sphere*sin(phi) + yc;   % This is the height.

fprintf('%15s  %24.16f\n','min(x)',min(x));
fprintf('%15s  %24.16f\n','max(x)',max(x));

cyl_x = [1,1,-1,-1,1]*R_cyl;
cyl_y = [-1,1,1,-1,-1]*H_cyl/2;

plot([-5,5],[0,0],'k--');
hold on;
plot([0,0],[-5,5],'k--');

plot(cyl_x,cyl_y,'k','linewidth',2);

th = linspace(0,2*pi,1000);
plot(R_sphere*cos(th) - L, R_sphere*sin(th), 'r', 'linewidth',2);
hold on;

plot(x,y,'b','linewidth',3);

% x(x < R_cyl) = nan;
% plot(x,y,'b','linewidth',4);

p(1) = plot(-L,0,'r.','markersize',35,'linewidth',30);
p(2) = plot(xbar,0,'g.','markersize',35,'linewidth',30);
lstr{1} = sprintf('L = %.4g\n',L);
lstr{2} = sprintf('xbar = %.4g\n',xbar);
legend(p,lstr,'fontsize',14);

axis([[-1,1]*R_cyl*4, [-1,1]*H_cyl/2*1.5]);
daspect([1,1,1]);

set(gcf,'position',[1649,386, 560, 420]);

end