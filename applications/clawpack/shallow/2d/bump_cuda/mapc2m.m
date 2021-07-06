function [xp,yp,zp] = mapc2m(xc,yc)

map = 'nomap';
% map = 'pillowdisk';

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'pillowdisk'
        [xp,yp,zp] = mapc2m_pillowdisk(xc,yc);
        xp = 2.5*xp;
        yp = 2.5*yp;
end
zp = 0*xp;

end