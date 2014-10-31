function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'brick';
% map = 'torus';
% map = 'latlong';

switch map
    case 'brick'
        mi = 4;
        ni = 4;
        isflat = true;        
        [xp,yp,zp] = mapc2m_brick(xc,yc);
        xp = xp/mi;
        yp = yp/ni;
        xp = 2*xp - 1;
        yp = 2*yp - 1;
    case 'torus'
        isflat = false;
        mi = 2;
        mj = 1;
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        xc1 = xc1/mi;
        yc1 = yc1/mj;
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
    case 'latlong'
        isflat = false;
        mi = 2;
        mj = 1;
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        xc1 = xc1/mi;
        yc1 = yc1/mj;
        a = 0.5-40/180;
        b = 0.5+40/180;
        yc1 = a + (b-a)*yc1;        
        [xp,yp,zp] = mapc2m_latlong(xc1,yc1);
end

end