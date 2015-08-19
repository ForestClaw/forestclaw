function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'cart';
% map = 'pillowdisk';
% map = 'squareddisk';
% map = 'pillowdisk5';
% map = 'fivepatch';

% This domain should be in [0,2],[0,2]

shift = [1,1,0];

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0.025;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
                
        xp = xp + shift(1);
        yp = yp + shift(2);
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        xp = xp + shift(1);
        yp = yp + shift(2);
        % (xp,yp) in [-1,1]x[-1,1]
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
        xp = xp + shift(1);
        yp = yp + shift(2);
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
        
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
        b = getblocknumber();
        s = 0.0;
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
        xp = xp+shift(1);
        yp = yp+shift(1);
end
zp = 0*xp;




end