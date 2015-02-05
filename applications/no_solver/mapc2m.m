function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'nomap';
% map = 'cart';
map = 'fivepatch';
map = 'pillowdisk';
map = 'pillowdisk5';

alpha = 0.5;

switch map
    case 'nomap'
        xp = xc;
        yp = yc;
        zp = 0*xc;
    case 'cart'
        [xp,yp,zp] = mapc2m_cart(xc,yc);
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
    case 'pillowdisk5'
        [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc);
end