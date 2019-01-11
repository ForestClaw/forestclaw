function [xp,yp,zp] = mapc2m_twisted_torus(xc,yc,alpha)

% alpha = 0.4;  % outer radius/inner radius

r = 1 + alpha*cos(2*pi*(xc+yc));
zp = alpha.*sin(2*pi*(xc+yc));
tw = 0;
th = 2*pi*(xc - tw*zp/0.4);
xp = r.*cos(th);
yp = r.*sin(th);

end