function [xp,yp,zp] = mapc2m(xc,yc)
global map_isflat;


map = 'torus';
% map = 'twisted_torus';
% map = 'flat';
% map = 'cart';


map_isflat = strcmp(map,'flat');
% R = 1;
% r = 0.4;

alpha = 0.4;

shift = [1,1,0];
scale = [0.5 0.5,1];

switch map
    case 'flat'
        s = 0;
        [xp,yp,~] = mapc2m_brick(xc,yc,s);
        zp = 0*xp;
         
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,~] = mapc2m_cart(xc1,yc1);
                
        xp = scale*xp;
        yp = scale*yp;
        xp = xp + shift(1);
        yp = yp + shift(2);
        zp = 0*xp;
        
    case 'torus'
        % Find center to expand radius around.
        s = 0.00;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1,alpha);

    case 'twisted_torus'
        % Find center to expand radius around.
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_twisted_torus(xc1,yc1,alpha);
end

end
