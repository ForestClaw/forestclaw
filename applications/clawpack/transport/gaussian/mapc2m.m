
function [xp,yp,zp] = mapc2m(xc1,yc1)

global notpillowsphere;

map = 'cubedsphere';
% map = 'pillowsphere';

switch map
    case 'nomap'
        % This needs work
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc1;
        yp = yc1;
        zp = zeros(size(xc1));
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
        %{
            case 2
                xp = xp + s;
            case 3
                yp = yp - s;
            case 4
                yp = yp + s;
            case 5
                zp = zp + s;
        end
        %}
end


end