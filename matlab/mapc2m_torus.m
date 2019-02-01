function [xp,yp,zp] = mapc2m_torus(xc,yc,alpha,beta)

if (nargin < 4)
    beta = 0;
end

r = alpha*(1 + beta*sin(2*pi*xc));
R = 1 + r.*cos(2*pi*yc);

% at yc = 0.5,  R = 1-r(x);
% X = (1-r(x))*cos(2*pi*x);
% Y = (1-r(x))*sin(2*pi*x);

xp = R.*cos(2*pi*xc);
yp = R.*sin(2*pi*xc);
zp = r.*sin(2*pi*yc);

end