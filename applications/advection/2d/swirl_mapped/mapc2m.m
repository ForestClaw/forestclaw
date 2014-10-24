function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'identity';
% map = 'cart';
% map = 'fivepatch';
map = 'brick';

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
    case 'brick'
        mi = 4;
        mj = 4;
        [xp,yp,zp] = mapc2m_brick(xc,yc);        
        xp = xp/mi;
        yp = yp/mj;
end
zp = 0*xp;


end