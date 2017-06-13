function [xp,yp,zp] = mapc2m(xc,yc)

global notpillowsphere;

map = 'pillowsphere5';
% map = 'pillow';

alpha = 0.4;

switch map
    case 'pillow'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);

    case 'pillowsphere5'
        notpillowsphere = true;
        [xc1,yc1,zc1] = mapc2m_fivepatch(0.5,0.5,alpha);
        xc1 = (xc1 + 1)/2;
        yc1 = (yc1 + 1)/2;
        [xp,yp,zp]= mapc2m_pillowsphere(xc1,yc1);
        v = [xp,yp,zp];

        [xp,yp,zp] = mapc2m_fivepatch(xc,yc,alpha);
        xc = (xp + 1)/2;
        yc = (yp + 1)/2;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
        s = 0.0;
        xp = xp + s*v(1);
        yp = yp + s*v(2);
        zp = zp + s*v(3);
end
