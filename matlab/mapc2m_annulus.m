function [xp,yp,zp] = mapc2m_annulus(xc,yc)

alpha = 0.4;
r = alpha + (1-alpha)*yc;
xp = r.*cos(2*pi*xc);
yp = r.*sin(2*pi*xc);
zp = 0*xp;

end