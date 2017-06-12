function [xp,yp,zp] = mapc2m(xc,yc)

map = 'fivepatch';
% map = 'cart';

switch map
    case 'cart'
        xp = xc;
        yp = yc;
        zp = 0*xp;
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
end


end