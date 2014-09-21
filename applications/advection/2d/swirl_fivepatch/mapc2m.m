function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'cart';
map = 'fivepatch';
% map = 'disk';
% map = 'pillow';

switch map
    case 'cart'
        xp = xc;
        yp = yc;
        zp = 0*xp;
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
    case 'disk'
        [xp1,yp1] = mapc2m_fivepatch(xc,yc);
        [xp,yp,zp] = mapc2m_pillowdisk(xp1,yp1);
        xp = (xp+1)/2;
        yp = (yp+1)/2;
    case 'pillow'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        xp = (xp+1)/2;
        yp = (yp+1)/2;
end


end