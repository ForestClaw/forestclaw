function [xp,yp,zp] = mapc2m(xc,yc)

beta = 0.4;
t1 = 0.125;
t2 = 0.375;
theta = [t1, t2];


[xp,yp,zp] = mapc2m_annulus(xc,yc,beta, theta);

end
