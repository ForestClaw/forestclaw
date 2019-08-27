function f = vel_ellipse(t,u)

global theta_range phi_range alpha_disk

pi2 = 2*pi;

[example,~,~,~, ~,~,~,~,theta_range,phi_range] = read_vars();

alpha = alpha_disk;

N = length(u)/3;

f = zeros(size(u));

X = u(1:N);
Y = u(N+1:2*N);
Z = u(2*N+1:3*N);

th0 = atan2(Y,X);
m = ~isreal(th0);
if sum(m(:)) > 0
    fprintf('vel_ellipse : th0 contains imaginary entries\n');
    keyboard;
end

m = th0 < 0;
th0(m) = th0(m) + pi2;

tr1 = pi2*theta_range(1);
tr2 = pi2*theta_range(2);
xm = (th0-tr1)/(tr2-tr1);

% m = Z/alpha > 1;
% Z(m) = alpha;

phi = asin(Z/alpha);
m = ~isreal(phi);
if sum(m(:)) > 0
    fprintf('vel_ellipse : phi contains imaginary entries\n');
    keyboard;
end

r = sqrt(X.^2 + Y.^2);

m1 = (r <= 1) & (phi > 0);
m2 = (r <= 1) & (phi <= 0);

phi(m1) =  pi - phi(m1);
phi(m2) = -pi - phi(m2);

m = phi < 0;
phi(m) = phi(m) + pi2;

% Assume phi in [0,2*pi]
ym = phi/pi2;

[t11,t21,t31, t12,t22,t32] = arrayfun(@covariant_basis,xm,ym);

if example == 0
    % Rigid body rotation
    u1 = -1;  % Revs per second
    u2 = 0;
    
elseif example == 1
    u1 = 0;
    u2 = 1;
end

vc1 = u1.*t11 + u2.*t12;
vc2 = u1.*t21 + u2.*t22;
vc3 = u1.*t31 + u2.*t32;

f(1:N) = vc1;
f(N+1:2*N) = vc2;
f(2*N+1:3*N) = vc3;
    
end

function [t11, t21, t31, t12, t22, t32] = covariant_basis(x,y)

global theta_range phi_range alpha_disk

pi2 = 2*pi;

alpha = alpha_disk;

tr1 = theta_range(1);
tr2 = theta_range(2);

pr1 = phi_range(1);
pr2 = phi_range(2);

theta = pi2*(tr1 +  (tr2-tr1)*x);
thetax = pi2*(tr2-tr1);

phi = pi2*(pr1 + (pr2-pr1)*y);
phiy = pi2*(pr2-pr1);

r = 1 + alpha*cos(phi);
ry = -alpha*sin(phi)*phiy;

t11 = -r.*sin(theta)*thetax;
t21 = r.*cos(theta)*thetax;
t31 = zeros(size(theta));
    
t12 = ry.*cos(theta);
t22 = ry.*sin(theta);
t32 = alpha*cos(phi)*phiy;

end
    
