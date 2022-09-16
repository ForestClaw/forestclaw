function [xp,yp,zp] = mapc2p(xc,yc,zc)

parms = read_vars();
if parms.example == 0
    map = 'nomap';
elseif parms.example == 1
    map = 'cart';   % brick
elseif parms.example == 2
    map = 'fivepatch';
elseif parms.example == 3
    map = 'bilinear';
end

% This domain should be in [0,1],[0,1]

alpha = parms.alpha;

scale = [0.5,0.5,1];
shift = [0.5,0.5,0];

switch map
    case 'nomap'
        xp = xc;
        yp = yc;
        zp = zc;
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,~] = mapc2m_cart(xc1,yc1);
                
        xp = scale(1)*xp + shift(1);
        yp = scale(2)*yp + shift(2);
        zp = zc;
        
    case 'fivepatch'
        [xp,yp,~] = mapc2m_fivepatch(xc,yc,alpha);

        xp = scale(1)*xp + shift(1);
        yp = scale(2)*yp + shift(2);
        zp = zc;

    case 'bilinear'
        center = parms.center;
        [xp,yp,~] = mapc2m_bilinear(xc,yc,center);

        xp = scale(1)*xp + shift(1);
        yp = scale(2)*yp + shift(2);
        zp = zc;

end




end