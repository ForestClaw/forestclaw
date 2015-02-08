function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

isflat = true;
s = 0.0;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
[xp,yp,zp] = mapc2m_cart(xc1,yc1);
xp = xp + 1;
yp = (yp + 1)/4;

zp = zeros(size(xp));

end