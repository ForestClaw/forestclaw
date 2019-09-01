function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,beta,theta,~,~] = read_vars();

s = 0;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
[xp,yp,zp] = mapc2m_annulus(xc1,yc1,beta, theta);

end
