function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'identity';
% map = 'cart';
map = 'fivepatch';


switch map
    case 'identity'
        xp = xc;
        yp = yc;        
    case 'cart'
        [xp,yp,zp] = mapc2m_cart(xc,yc);
        xp = (xp+1)/2;
        yp = (yp+1)/2;
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
        xp = (xp+1)/2;
        yp = (yp+1)/2;
end
zp = 0*xp;


end