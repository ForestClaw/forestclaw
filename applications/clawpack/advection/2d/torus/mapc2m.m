function [xp,yp,zp] = mapc2m(xc,yc)
global map_isflat;


map = 'torus';
map = 'twisted_torus2';
% map = 'flat';
% map = 'cart';

switch map
    case 'torus'
       L = eye(2);
    case 'twisted_torus'
       L = [1 0; 1 1];
    case 'twisted_torus2'
       L = [1 1; 1 0];
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;


map_isflat = strcmp(map,'flat');
% R = 1;
% r = 0.4;

alpha = 0.4;
beta = 0;

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
        
    case {'torus','twisted_torus','twisted_torus2'}
        % Find center to expand radius around.
        s = 0.00;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        
        [xp,yp,zp] = mapc2m_torus(a(xc1,yc1),b(xc1,yc1),alpha,beta);
        % [xp,yp,zp] = mapc2m_torus(xc1,yc1,alpha,beta);

end

end
