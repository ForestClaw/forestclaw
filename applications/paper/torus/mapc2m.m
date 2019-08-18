function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,alpha,beta,~,~,tr,pr] = read_vars();


% Find center to expand radius around.
s = 0.00;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
xc2 = tr(1) + (tr(2)-tr(1))*xc1;
yc2 = pr(1) + (pr(2)-pr(1))*yc1;
[xp,yp,zp] = mapc2m_torus(xc2,yc2,alpha,beta);

end
