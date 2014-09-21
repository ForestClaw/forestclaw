function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'squareddisk';
% map = 'pillowdisk';
map = 'cart';
% map = 'cart_nomap';

switch map
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        xp = xp + 1;
        yp = yp + 1;
    case 'squareddisk'
        [xp,yp,zp] = mapc2m_squareddisk(xc,yc);
        xp = xp + 1;
        yp = yp + 1;        
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        xp = 2*xc;
        yp = 2*yc;
        zp = 0*xp;
    case 'cart_nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        xp = xc;
        yp = yc;
        zp = 0*xp;
end



end