function [xp,yp,zp] = mapc2m_torus(xc,yc,alpha,beta)

if (nargin < 4)
    beta = 0;
end

pi2 = 2*pi;

r = alpha*(1 + beta*sin(pi2*xc));
R = 1 + r.*cos(pi2*yc);

xp = R.*cos(pi2*xc);
yp = R.*sin(pi2*xc);
zp = r.*sin(pi2*yc);

end