function [xp,yp,zp] = mapc2m(xc,yc)

[~,~,~,R, H,~,~,~,] = read_vars();

[xp,yp,zp] = mapc2m_cylinder(xc,yc,R,H);
   
end
