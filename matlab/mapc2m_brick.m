function [xp,yp,zp] = mapc2m_brick(xc,yc)

% (mi,mj) brick

ix = [0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3];
jy = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3];

blockno = getblocknumber();

j = jy(blockno+1);
i = ix(blockno+1);
xp = i + xc;
yp = j + yc;
zp = 0*xp;

end