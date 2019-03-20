function [xp,yp,zp] = mapc2m_annulus(xc,yc,beta)

r = beta + (1-beta)*yc;
xp = r.*cos(2*pi*xc);
yp = r.*sin(2*pi*xc);
zp = 0*xp;

end