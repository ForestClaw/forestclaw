
function [xp,yp,zp] = mapc2m(xc1,yc1)

global notpillowsphere flat;

flat = false;

map_list = {'cubedsphere', 'pillowsphere'};

[example,mapping,ic,~, ~] = read_vars();

map = map_list{mapping+1};

% map = 'cubedsphere_flat';

switch map
    case 'pillowsphere'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc1,yc1);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc1,yc1);
        b = getblocknumber();
        s = 0;
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
    case 'cubedsphere_flat'
        flat = true;
        [xp,yp,zp] = mapc2m_cubedsphere_flat(xc1,yc1);
end


end