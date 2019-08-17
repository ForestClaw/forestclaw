function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,alpha,beta,~] = read_vars();


% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
[xp,yp,zp] = mapc2m_torus(xc1,yc1,alpha,beta);

end
