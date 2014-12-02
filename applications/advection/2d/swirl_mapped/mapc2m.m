function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'identity';
% map = 'cart';
map = 'fivepatch';

switch map
    case 'identity'
        xp = xc;
        yp = yc;        
    case 'cart'
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
        
        brick_data = load('brick.dat');
        mi = brick_data(1,1);
        mj = brick_data(1,2);
        xv = brick_data(2:end,1);
        yv = brick_data(2:end,2);

        blockno = getblocknumber();
        s = 0.02;
        xs = (xv(blockno+1)+1-mi/2)*s;
        ys = (yv(blockno+1)+1-mj/2)*s;
        xp = xp + xs;
        yp = yp + ys;

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
        % xp = (xp+1)/2;
        % yp = (yp+1)/2;
        
    case 'brick'
        [xp,yp,zp] = mapc2m_brick(xc,yc);        
end
zp = 0*xp;


end