function [xp,yp,zp] = mapc2m(xc,yc)

d = load('_output/mapping.dat');

beta = d(5);
theta = d(6:7);


[xp,yp,zp] = mapc2m_annulus(xc,yc,beta, theta);

end
