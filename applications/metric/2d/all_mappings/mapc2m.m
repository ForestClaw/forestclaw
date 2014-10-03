function [xp,yp,zp] = mapc2m(xc,yc)

global notpillowsphere;


% map = 'cart';            % 1
% map = 'fivepatch';       % 2
% map = 'squareddisk';     % 3
% map = 'pillowdisk';      % 4
% map = 'pillowdisk5';     % 5
map = 'pillowsphere';    % 6
% map = 'cubedsphere';     % 7
% map = 'torus';           % 8

alpha = 0.5;

notpillowsphere = true;
switch map
    case 'nomap'
        xp = xc;
        yp = yc;
    case 'cart'
        [xp,yp,zp] = mapc2m_cart(xc,yc);
    case 'fivepatch'
        [xp,yp,zp] = mapc2m_fivepatch(xc,yc);
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
    case 'squareddisk'
        [xp,yp,zp] = mapc2m_squareddisk(xc,yc);
    case 'pillowdisk5'
        [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc);
    case 'pillowsphere'
        notpillowsphere = false;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc,yc);
    case 'torus'
        [xp,yp,zp] = mapc2m_torus(xc,yc);
        
end
       
        
end