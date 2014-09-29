function [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc)

[xc,yc,zc] = mapc2m_fivepatch(xc,yc);

xc = (xc+1)/2;
yc = (yc + 1)/2;

[xp,yp,zp] = mapc2m_pillowdisk(xc,yc);

end