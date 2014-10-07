function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'cart';
% map = 'fivepatch';

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
end