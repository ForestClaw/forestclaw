function [xp,yp,zp] = mapc2m(xc,yc)

% map = 'nomap';
% map = 'cart';
% map = 'fivepatch';
% map = 'pillowdisk';
% map = 'squareddisk';
% map = 'pillowdisk5';
map = 'pillowsphere';
% map = 'cubedsphere';


alpha = 0.5;

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
        [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc);
    case 'pillowsphere'
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc,yc);
end
       
        
end