
function [xp,yp,zp] = mapc2m(xc1,yc1)

global notpillowsphere;

map = 'cubedsphere';
% map = 'pillowsphere';

b = getblocknumber();

scale = [1.5, 1.5, 1.5];

switch map
    case 'pillowsphere'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc1,yc1);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc1,yc1);
        s = 0.0;
        switch b
            case 0
                zp = zp - s;
            case 1
                yp = yp + s;
            case 2
                xp = xp - s;
            case 3
                zp = zp + s;
            case 4
                yp = yp - s;
            case 5
                xp = xp + s;
            otherwise
        end
        xp = scale(1)*xp;
        yp = scale(2)*yp;
        zp = scale(3)*zp;
end


end