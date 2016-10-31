function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'pillowdisk5';

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'pillowdisk5'
        [xp,yp,zp] = mapc2m_pillowdisk5(xc,yc);
end
zp = 0*xp;

end
