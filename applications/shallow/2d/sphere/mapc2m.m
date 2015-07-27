function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

% map = 'latlong';
map = 'cubedsphere';

isflat = false;

R = 1;
r = 0.4;

ll = load('latlong.dat');
lat = ll(1:2);
lng = ll(3:4);



switch map
    case 'latlong'
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);

        % Map into [0,1]x[0,1]
        xc2 = lng(1) + (lng(2) - lng(1))*xc1;
        yc2 = lat(1) + (lat(2) - lat(1))*yc1;
        [xp,yp,zp] = mapc2m_latlong(xc2,yc2);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc,yc);

end

if (isflat)
    zp = zeros(size(xp));
end

end
