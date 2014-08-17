function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'squareddisk';
% map = 'pillowdisk';
map = 'cart';

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
        xp = 2*xc;
        yp = 2*yc;
        zp = 0*xp;
end



end