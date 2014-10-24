function [xp,yp,zp] = mapc2m_brick(xc,yc)

% (mi,mj) brick

blockno = getblocknumber();

mi = 4;
mj = 4;

j = floor(blockno/mi);
i = blockno - j*mi;
xp = i + xc;
yp = j + yc;
zp = 0*xp;

end