
function [xp,yp,zp] = mapc2m(xc1,yc1)

global notpillowsphere;

% map = 'pillowsphere';
map = 'cubedsphere';

switch map
    case 'pillowsphere'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc1,yc1);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc1,yc1);
end


end