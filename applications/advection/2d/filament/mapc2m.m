function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'nomap';
map = 'cart';
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
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        % (xp,yp) in [-1,1]x[-1,1]
    case 'squareddisk'
        [xp,yp,zp] = mapc2m_squareddisk(xc,yc);
        % (xp,yp) in [-1,1]x[-1,1]
    case 'pillowdisk5'
        [xp1,yp1,zp] = mapc2m_fivepatch(xc,yc);
        xp1 = (xp1 + 1)/2;
        yp1 = (yp1 + 1)/2;
        [xp,yp,zp] = mapc2m_pillowdisk(xp1,yp1);
end
if (strcmp(map,'nomap') == 0)
    xp = xp + shift(1);
    yp = yp + shift(2);
end
zp = 0*xp;




end