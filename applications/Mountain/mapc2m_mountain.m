function [xp,yp,zp] = mapc2m_mountain(xc,yc,s,x_scale,y_scale)

% Shrink block
xc = (1-s)*xc + s/2;
yc = (1-s)*yc + s/2;

brick_data = load('brick.dat');
mi = brick_data(1,1);
mj = brick_data(1,2);
xv = brick_data(2:end,1);
yv = brick_data(2:end,2);

blockno = getblocknumber();

% get xp
xp = x_scale*(xv(blockno+1) + xc)/mi;

% Get yp
xl = xv(blockno+1);
yl = yv(blockno+1);

% Map corners to physical space
xp_corner = x_scale*[xl xl+1 xl+1 xl]/mi;
yp_corner = y_scale*[yl yl yl+1 yl+1]/mj;

height_corners = mountain_height(xp_corner);

f = yp_corner/y_scale;
yp_corner = height_corners + f.*(y_scale-height_corners);

ztop = yp_corner(4) + xc*(yp_corner(3) - yp_corner(4));
if (yl == 0)
    zbot = mountain_height(xp);
else
    zbot = yp_corner(1) + xc*(yp_corner(2) - yp_corner(1));
end
yp = zbot + yc.*(ztop-zbot);


zp = 0*xp;



% Get four corners of physical block

% decide if block is at the bottom edge or not

% Map (xc,yc) \in [0,1]x[0,1] into block


end