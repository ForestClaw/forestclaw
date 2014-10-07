function [xp,yp,zp] = mapc2m(xc,yc)

global notpillowsphere;

map = 'pillow';
% map = 'pillowsphere5';

switch map
    case 'pillow'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
    case 'pillowsphere5'
        notpillowsphere = true;
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);  
        xc = (xp + 1)/2;
        yc = (yp + 1)/2;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
end