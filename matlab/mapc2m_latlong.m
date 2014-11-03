function [xp,yp,zp] = mapc2m(xc,yc)

% assume that [xc,yc] in [0,mi], [0,mj]

xp = cos(pi*(0.5-yc)).*cos(2*pi*xc);
yp = cos(pi*(0.5-yc)).*sin(2*pi*xc);
zp = sin(pi*(0.5-yc));

end


