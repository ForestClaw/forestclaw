function [xp,yp,zp] = mapc2m(xc,yc)

% ---------------------------------------------------------
% Assume a unit box (xc,yc) \in [-1,1]x[-1,1]
% ---------------------------------------------------------

[xp,yp] = mapc2m_inclusions(xc,yc);

zp = 0*xp;
