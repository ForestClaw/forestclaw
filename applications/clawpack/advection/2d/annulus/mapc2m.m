function [xp,yp,zp] = mapc2m(xc,yc)

global square;

map = 'twisted_annulus';

square = false;
if (square)
    map = 'cart';
end

twist = load('twist.dat');
beta = 0.4;

switch map
    case {'annulus','cart'}
       L = eye(2);
    case 'twisted_annulus'
       L = [1 twist; 0 1];
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;


s = 0;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);

switch map
    case 'cart'
        square = true;
        xp = a(xc1,yc1);
        yp = b(xc1,yc1);
        zp = 0*xc1;
    otherwise
        [xp,yp,zp] = mapc2m_annulus(a(xc1,yc1),b(xc1,yc1),beta);
end

% isflat = true;
% s = 0.0;
% [xp,yp,zp] = mapc2m_annulus(a(xc1,yc1),b(xc1,yc1),beta);
% [xp,yp,zp] = mapc2m_mobius(a(xc1,yc1),b(xc1,yc1));

end
