function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'cart';
% map = 'torus';
map = 'latlong';

switch map
    case 'cart'
        isflat = true;        
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
    case 'torus'
        isflat = false;
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
    case 'latlong'
        isflat = false;
        lat = [-70 70];
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        a = 0.5+lat(1)/180;
        b = 0.5+lat(2)/180;
        yc1 = a + (b-a)*yc1;        
        [xp,yp,zp] = mapc2m_latlong(xc1,yc1);
end

end