function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

% map = 'nomap';
% map = 'brick';
map = 'torus';

switch map
    case 'nomap'
        isflat = true;
        xp = xc;
        yp = yc;
        zp = 0*xp;
    case 'brick'
        mi = 2;
        ni = 2;
        isflat = true;        
        [xp,yp,zp] = mapc2m_brick(xc,yc);
        xp = xp - 1;
        yp = yp - 1;
    case 'torus'
        isflat = false;
        mi = 2;
        mj = 1;
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        xc1 = xc1/mi;
        yc1 = yc1/mj;
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
end

end