function [xp,yp,zp] = mapc2m(xc,yc)

map = 'identity';
% map = 'cart';
% map = 'fivepatch';

switch map
    case 'identity'
        xp = xc;
        yp = yc;        
    case 'cart'
        s = 0.05;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
        
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
        b = getblocknumber();
        s = 0.025;
        switch b
            case 0
                yp = yp - s;
            case 1
                xp = xp - s;
            case 3
                xp = xp + s;
            case 4
                yp = yp + s;
        end
        xp = (xp+1)/2;
        yp = (yp+1)/2;
        
end
zp = 0*xp;


end