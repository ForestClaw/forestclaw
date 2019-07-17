function [xp,yp,zp] = mapc2m(xc,yc)

global square;

% map = 'twisted_annulus';
map = 'annulus';

square = false;
if (square)
    map = 'cart';
end

twist = load('twist.dat');
beta = 0.4;

t1 = 0.125;
t2 = 0.375;
theta = [t1, t2];


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
        [xp,yp,zp] = mapc2m_annulus(a(xc1,yc1),b(xc1,yc1),beta,theta);
end

end
