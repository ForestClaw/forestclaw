function [xp,yp,zp] = mapc2p(xc,yc,zc)

parms = read_vars();

map = 'latlong';

switch map
    case 'latlong'
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);

        % Map into [0,1]x[0,1]
        lat = parms.latitude;
        lng = parms.longitude;
        xc2 = lng(1) + (lng(2) - lng(1))*xc1;
        yc2 = lat(1) + (lat(2) - lat(1))*yc1;
        [xp1,yp1,zp1] = mapc2m_latlong(xc2,yc2);

        phi = asin(zp1);          % returns value in [-pi/2, pi/2]
        theta = atan2(yp1,xp1);    % returns value in [-pi, pi]
        m = theta < 0;
        theta(m) = theta(m) + 2*pi;

        % Assume zc in [0,1]
        R = parms.maxelev*zc + 1;    
        xp = R.*cos(phi).*cos(theta);
        yp = R.*cos(phi).*sin(theta);
        zp = R.*sin(phi);
        
end

end
