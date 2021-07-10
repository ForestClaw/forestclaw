function [xp,yp,zp] = mapc2m(xc,yc)

d = load('setprob.data');
alpha = d(2);
beta = d(3);

% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
[xp,yp,zp] = mapc2m_torus(xc1,yc1,alpha, beta);


end
