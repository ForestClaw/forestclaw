function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

% map = 'nomap';
% map = 'cart';
map = 'torus';

switch map
    case 'nomap'
        isflat = true;
        xp = xc;
        yp = yc;
        zp = 0*xp;
    case 'cart'
        isflat = true;
        [xp,yp,zp] = mapc2m_cart(xc,yc);
    case 'torus'
        isflat = false;
        [xp,yp,zp] = mapc2m_torus(xc,yc);
end

end