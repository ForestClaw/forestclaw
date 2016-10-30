function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'latlong';

R = 1;
r = 0.4;


switch map
    case 'nomap'
        isflat = true;
        xp = xc;
        yp = yc;

    case 'brick'
        isflat = true;
        s = 0.0;
        [xp,yp,~] = mapc2m_brick(xc,yc,s);
         b = load('brick.dat');
%          mi = b(1,1);
%          mj = b(1,2);
%          xp = mi*xp;
%          yp = mj*yp;

    case 'cart'
        isflat = true;
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);

    case 'torus'
        isflat = false;
        % Find center to expand radius around.
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
        
    case 'latlong'
        isflat = false;
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);

        % Map into [0,1]x[0,1]
        lng = [0 360];
        lat = [-50 50];
        xc2 = lng(1) + (lng(2) - lng(1))*xc1;
        yc2 = lat(1) + (lat(2) - lat(1))*yc1;
        [xp,yp,zp] = mapc2m_latlong(xc2,yc2);
    case 'annulus'
        isflat = true;
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_annulus(xc1,yc1);
    case 'duplicate'
        isflat = true;
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
%         b = load('brick.dat');
%         mi = b(1,1);
%         mj = b(1,2);
        xp = -1 + 2*xc1;
        yp = -1 + 2*yc1;
end

if (isflat)
    zp = zeros(size(xp));
end

end
