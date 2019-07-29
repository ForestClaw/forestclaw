function [xp,yp,zp] = mapc2m(xc,yc)

global square

map = 'twisted_annulus';
% map = 'cart';
square = true;

beta = 0.4;
twist = load('twist.dat');

switch map
    case {'cart','annulus'}
       L = eye(2);
    case 'twisted_annulus'
       L = [1 twist; 0 1];
    otherwise
        % Nothing
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;

square = false;
switch map
    case 'cart'
        square = true;
        xp = a(xc,yc);
        yp = b(xc,yc);
        zp = 0*xc;
    otherwise
        [xp,yp,zp] = mapc2m_annulus(a(xc,yc),b(xc,yc),beta);
end

end
