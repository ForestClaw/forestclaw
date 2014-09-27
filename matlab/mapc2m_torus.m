function [xp,yp,zp] = mapc2m_torus(xc,yc)

alpha = 0.4;

d = 1 + alpha*cos(2*pi*yc);
xp = d.*cos(2*pi*xc);
yp = d.*sin(2*pi*xc);
zp = alpha.*sin(2*pi*yc);

end