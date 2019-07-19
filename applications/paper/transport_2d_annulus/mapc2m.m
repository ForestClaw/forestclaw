function [xp,yp,zp] = mapc2m(xc,yc)

d = load('_output/mapping.dat');

beta = d(3);
t1 = d(4);
t2 = d(5);
theta = [t1, t2];


[xp,yp,zp] = mapc2m_annulus(xc,yc,beta, theta);

end
