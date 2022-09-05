function [xp,yp,zp] = mapc2p(xc,yc,zc)

global notpillowsphere

parms = read_vars();

map_list = {'nomap', 'cart', 'latlong', 'cubedsphere', 'pillowsphere'};

map = map_list{parms.mapping+1};

scale = [0.5,0.5,0.5];
shift = [0.5,0.5,0.5];

switch map
    case 'nomap'
        xp = xc;
        yp = yc;
        zp = zc;
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,~] = mapc2m_cart(xc1,yc1);                
        zp = -1 + 2*zc;

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

    case {'cubedsphere','pillowsphere'}
        if strcmp('pillowsphere',map) == 1
            notpillowsphere = false;
            [xp1,yp1,zp1] = mapc2m_pillowsphere(xc,yc);
        else
            [xp1,yp1,zp1] = mapc2m_cubedsphere(xc,yc);
        end

        phi = asin(zp1);          % returns value in [-pi/2, pi/2]
        theta = atan2(yp1,xp1);    % returns value in [-pi, pi]
        m = theta < 0;
        theta(m) = theta(m) + 2*pi;

        % Assume zc1 in [0,1]
        R = parms.maxelev*zc + 1;
        xp = R.*cos(phi).*cos(theta);
        yp = R.*cos(phi).*sin(theta);
        zp = R.*sin(phi);

end

end
