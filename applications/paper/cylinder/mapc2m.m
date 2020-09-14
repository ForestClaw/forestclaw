function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,R,H,~,~,~,tr,pr] = read_vars();


% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
th = 2*pi*xc1;
xp = R*cos(th);
yp = R*sin(th);
zp = H*yc1;

% [xp,yp,zp] = mapc2m_torus(xc2,yc2,alpha,beta);

end
