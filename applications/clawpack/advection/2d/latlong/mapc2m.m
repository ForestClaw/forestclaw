function [xp,yp,zp] = mapc2m(xc,yc)

map = 'latlong';

scale = [1.5, 1.5, 1.5];

switch map
    case 'latlong'
        isflat = false;
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);

        % Map into [0,1]x[0,1]
        % Map into [0,1]x[0,1]
        parms = read_vars();
        lat = parms.latitude;
        lng = parms.longitude;
        xc2 = lng(1) + (lng(2) - lng(1))*xc1;
        yc2 = lat(1) + (lat(2) - lat(1))*yc1;
        [xp,yp,zp] = mapc2m_latlong(xc2,yc2);
        xp = scale(1)*xp;
        yp = scale(2)*yp;
        zp = scale(3)*zp;
end

end
