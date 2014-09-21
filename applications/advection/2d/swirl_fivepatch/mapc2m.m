function [xp,yp,zp] = mapc2m(xc,yc)

map = 'fivepatch';
% map = 'nomap';

switch map
    case 'nomap'
        xp = xc;
        yp = yc;
        zp = 0*xp;
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
end


end