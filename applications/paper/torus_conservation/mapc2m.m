function [xp,yp,zp] = mapc2m(xc,yc)

alpha = 0.4;
beta = 0.6;

% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);

% [xp,yp,zp] = mapc2m_torus(a(xc1,yc1),b(xc1,yc1),alpha,beta);
[xp,yp,zp] = mapc2m_torus(xc1,yc1,alpha,beta);

end
