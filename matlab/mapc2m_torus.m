function [xp,yp,zp] = mapc2m_torus(xc,yc,alpha)

r = 1 + alpha*cos(2*pi*yc);
zp = alpha.*sin(2*pi*yc);
th = 2*pi*xc;
xp = r.*cos(th);
yp = r.*sin(th);

end