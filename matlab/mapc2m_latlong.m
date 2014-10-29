function [xp,yp,zp] = mapc2m(xc,yc)

% assume that [xc,yc] in [0,mi], [0,mj]

xp = cos(pi*(yc-0.5)).*cos(2*pi*xc);
yp = cos(pi*(yc-0.5)).*sin(2*pi*xc);
zp = sin(pi*(yc-0.5));

end


