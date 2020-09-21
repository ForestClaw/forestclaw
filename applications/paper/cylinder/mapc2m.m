function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,R,H] = read_vars();


% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);

theta = 2*pi*xc1;
z = yc1;


xp = R*cos(theta);
yp = R*sin(theta);
zp = H*z;

% xp = R.*cos(theta);
% zp = R.*sin(theta);
% yp = H*z;
% zp(zp < 0) = nan;

end
