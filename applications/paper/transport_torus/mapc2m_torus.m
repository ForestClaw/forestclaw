function [xp,yp,zp] = mapc2m_torus2(xc,yc,alpha,beta,tr,pr)

if (nargin < 4)
    beta = 0;
end

pi2 = 2*pi;

theta = pi2*(tr(1) + (tr(2)-tr(1))*xc);
phi = pi2*(pr(1) + (pr(2) - pr(1))*yc);

r = alpha*(1 + beta*sin(theta));
R = 1 + r.*cos(phi);

xp = R.*cos(theta);
yp = R.*sin(theta);
zp = r.*sin(phi);

end