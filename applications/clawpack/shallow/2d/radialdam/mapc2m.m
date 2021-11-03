function [xp,yp,zp] = mapc2m(xc,yc)

% Choice 0 : nomap
% Choice 1 : five patch (square)
% Choice 2 : pillowdisk
% Choice 3 : pillowdisk5 (five patch square --> disk)

map_choice = 1;
map_list = {'nomap', 'fivepatch','pillowdisk','pillowdisk5'};

map = map_list{map_choice + 1};
% map = 'nomap';
% map = 'pillowdisk';
% map = 'pillowdisk5';
% map = 'fivepatch';

switch map
    case 'nomap'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        xp = xc;
        yp = yc;
    case 'fivepatch'
        alpha = 0.4;
        [xp,yp,~] = mapc2m_fivepatch(xc,yc,alpha);
        xp = 2.5*xp;
        yp = 2.5*yp;
    case 'pillowdisk'
        [xp,yp,~] = mapc2m_pillowdisk(xc,yc);
        xp = 2.5*xp;
        yp = 2.5*yp;
    case 'pillowdisk5'
        alpha = 0.4;
        [xp,yp,~] = mapc2m_pillowdisk5(xc,yc,alpha);
        xp = 2.5*xp;
        yp = 2.5*yp;

end
zp = 0*xp;

end