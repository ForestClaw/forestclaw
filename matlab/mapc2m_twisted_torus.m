function [xp,yp,zp] = mapc2m_twisted_torus(xc,yc)

alpha = 0.4;  % outer radius/inner radius

r = 1 + alpha*cos(2*pi*(xc+yc));
xp = r.*cos(2*pi*xc);
yp = r.*sin(2*pi*xc);
zp = alpha.*sin(2*pi*(xc+yc));

end