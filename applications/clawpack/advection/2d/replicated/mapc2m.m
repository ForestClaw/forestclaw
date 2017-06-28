function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'nomap';

R = 1;
r = 0.4;


switch map
    case 'nomap'
        isflat = true;
        xp = xc;
        yp = yc;

    case 'brick'
        isflat = true;
        s = 0.01;
        [xp,yp,~] = mapc2m_brick(xc,yc,s);
        b = load('brick.dat');
%         mi = b(1,1);
%         mj = b(1,2);
%         xp = mi*xp;
%         yp = mj*yp;
end

zp = zeros(size(xp));

end
