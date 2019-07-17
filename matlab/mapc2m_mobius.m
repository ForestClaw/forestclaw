function [xp,yp,zp] = mapc2m_mobius(x,y)

% Convert x in [0,1] to u in [0,2*pi]
% Convert y in [0,1] to v in [-1,1]

u = 2*pi*x;
v = -1 + 2*y;

r1 = v/2;
r = 1 + r1.*cos(u/2);
xp = r.*cos(u);
yp = r.*sin(u);
zp = r1.*sin(u/2);


end