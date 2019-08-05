function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,~,beta,theta,~,~] = read_vars();

[xp,yp,zp] = mapc2m_annulus(xc,yc,beta, theta);

end
