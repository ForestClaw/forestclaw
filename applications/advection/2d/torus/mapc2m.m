function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'torus';
% map = 'twisted_torus';


R = 1;
r = 0.4;


switch map
    case 'nomap'
        isflat = true;
        xp = xc;
        yp = yc;

    case 'brick'
        isflat = true;
        s = 0.05;
        [xp,yp,~] = mapc2m_brick(xc,yc,s);
         b = load('brick.dat');
%          mi = b(1,1);
%          mj = b(1,2);
%          xp = mi*xp;
%          yp = mj*yp;

    case 'torus'
        isflat = false;
        % Find center to expand radius around.
        s = 0.00;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);

    case 'twisted_torus'
        isflat = false;
        % Find center to expand radius around.
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_twisted_torus(xc1,yc1);


end

if (isflat)
    zp = zeros(size(xp));
end

end
