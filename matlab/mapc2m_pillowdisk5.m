function [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc)

% Returns map in [-1,1]x[-1,1]
[xc,yc,zc] = mapc2m_fivepatch(xc,yc);

% Map to [0,1]x[0,1]
xc = (xc + 1)/2;
yc = (yc + 1)/2;

[xp,yp,zp] = mapc2m_pillowdisk(xc,yc);

end