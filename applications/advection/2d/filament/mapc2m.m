function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'cart';
% map = 'pillowdisk';
% map = 'squareddisk';
% map = 'pillowfivepatch';

shift = [1,1,0];

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        [xp,yp,zp] = mapc2m_cart(xc,yc);
        xp = xp + shift(1);
        yp = yp + shift(2);
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        % (xp,yp) in [-1,1]x[-1,1]
    case 'squareddisk'
        [xp,yp,zp] = mapc2m_squareddisk(xc,yc);
        s = 0.05;
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
    case 'pillowfivepatch'
        [xp1,yp1,zp] = mapc2m_fivepatch(xc,yc);
        xp1 = (xp1 + 1)/2;
        yp1 = (yp1 + 1)/2;
        [xp,yp,zp] = mapc2m_pillowdisk(xp1,yp1);
        s = 0.05;
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