function [xp,yp,zp] = mapc2m_annulus(xc,yc,beta,theta)

t1 = theta(1);
t2 = theta(2);

r = beta + (1-beta)*yc;
t = t1 + (t2-t1)*xc;
xp = r.*cos(2*pi*t);
yp = r.*sin(2*pi*t);
zp = 0*xp;

end