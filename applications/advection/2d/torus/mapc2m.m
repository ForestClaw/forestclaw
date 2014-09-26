function [xp,yp,zp] = mapc2m(xc,yc)

global iscart;

map = 'cart';
% map = 'torus';

switch map
    case 'cart'
        iscart = true;
        [xp,yp,zp] = mapc2m_cart(xc,yc);
    case 'torus'
        iscart = false;
        [xp,yp,zp] = mapc2m_torus(xc,yc);
end

end