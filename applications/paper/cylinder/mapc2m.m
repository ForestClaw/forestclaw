function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,R_cyl,H_cyl] = read_vars();


% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);



theta = 2*pi*xc1;
z = yc1;


xp = R_cyl*cos(theta);
yp = R_cyl*sin(theta);
zp = H_cyl*z;

% xp = R.*cos(theta);
% zp = R.*sin(theta);
% yp = H*z;
% zp(zp < 0) = nan;

end
