function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

% map = 'cart';
map = 'torus';
% map = 'latlong';
% map = 'annulus';

R = 1;
r = 0.4;

switch map
    case 'cart'
        isflat = true;        
        s = 0.05;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
        
    case 'torus'
        isflat = false;
        % Find center to expand radius around.    
        s = 0.05;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
    case 'latlong'
        isflat = false;
        lat = [-50 50];
        a = 0.5+lat(1)/180;
        b = 0.5+lat(2)/180;
        s = 0.05;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        yc1 = a + (b-a)*yc1;        
        [xp,yp,zp] = mapc2m_latlong(xc1,yc1);
    case 'annulus'
        isflat = true;
        s = 0.05;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_annulus(xc1,yc1);
end

end