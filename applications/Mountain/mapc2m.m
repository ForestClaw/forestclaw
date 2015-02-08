function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'cutcell';
map = 'terrain';

x_scale = 5000;
y_scale = 2000;

s = 0.0;
switch map
    case 'cutcell'
        [xp,yp,zp] = mapc2m_brick(xc,yc,s);
        xp = x_scale*xp;
        yp = y_scale*yp;
    case 'terrain' 
        s = 0.05;
        [xp,yp,zp] = mapc2m_mountain(xc,yc,s,x_scale,y_scale);
end
zp = 0*xp;
end