function [xp,yp,zp] = mapc2m(xc,yc)

global notpillowsphere;

map = 'pillow';
% map = 'fivepatch';

switch map
    case 'pillow'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
    case 'fivepatch'
        notpillowsphere = true;
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);  
        [xp,yp,zp] = mapc2m_pillowsphere(xp,yp);
end