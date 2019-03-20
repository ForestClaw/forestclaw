function [xp,yp,zp] = mapc2m(xc,yc)

map = 'twisted_annulus';

switch map
    case 'annulus'
       L = eye(2);
    case 'twisted_annulus'
       L = [1 -0.2; 0 1];
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;


beta = 0.2;

isflat = true;
s = 0.0;
[xc1,yc1,~] = mapc2m_brick(xc,yc,s);
[xp,yp,zp] = mapc2m_annulus(a(xc1,yc1),b(xc1,yc1),beta);
% [xp,yp,zp] = mapc2m_mobius(a(xc1,yc1),b(xc1,yc1));

%xp = a(xc1,yc1);
%yp = b(xc1,yc1);

end
