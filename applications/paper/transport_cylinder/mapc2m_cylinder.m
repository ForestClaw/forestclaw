function [xp,yp,zp] = mapc2m_cylinder(xc,yc,R,H)


[~,~,~, R, H,~] = read_vars();

pi2 = 2*pi;

theta = pi2*xc;
z = yc;

xp = R.*cos(theta);
yp = R.*sin(theta);
zp = H*z;

% xp = R.*cos(theta);
% zp = R.*sin(theta);
% yp = H*z;
% zp(zp < 0) = nan;



end