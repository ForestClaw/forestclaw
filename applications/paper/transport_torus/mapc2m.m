function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,alpha, beta,~, ...
    ~,~,theta_range,phi_range] = read_vars();

[xp,yp,zp] = mapc2m_torus(xc,yc,alpha,beta,theta_range,phi_range);
   
end
