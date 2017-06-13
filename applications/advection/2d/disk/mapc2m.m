function [xp,yp,zp] = mapc2m(xc,yc)

map = 'pillowdisk';
map = 'pillowdisk5';

% This domain should be in [0,2],[0,2]

shift = [1,1,0];

switch map
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        xp = xp + shift(1);
        yp = yp + shift(2);
        % (xp,yp) in [-1,1]x[-1,1]
    case 'pillowdisk5'
        [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc);
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
        xp = xp + shift(1);
        yp = yp + shift(2);
        
end
zp = 0*xp;




end