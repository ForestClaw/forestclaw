function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'torus';
map = 'twisted_torus';
% map = 'flat';


isflat = false;
R = 1;
r = 0.4;


switch map
    case 'flat'
        s = 0;
        [xp,yp,~] = mapc2m_brick(xc,yc,s);
        zp = 0*xp;
         
    case 'torus'
        % Find center to expand radius around.
        s = 0.00;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);

    case 'twisted_torus'
        % Find center to expand radius around.
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,zp] = mapc2m_twisted_torus(xc1,yc1);
end

end
