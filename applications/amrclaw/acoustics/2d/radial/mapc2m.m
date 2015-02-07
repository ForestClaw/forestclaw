function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'squareddisk';

% This domain should be in [0,2],[0,2]

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'squareddisk'
        [xp,yp,zp] = mapc2m_squareddisk(xc,yc);
        s = 0.0;
        b = getblocknumber();
        switch b
            case 0
                yp = yp - s;
            case 1
                xp = xp - s;
            case 3
                xp = xp + s;
            case 4
                yp = yp + s;
        end
end
zp = 0*xp;




end