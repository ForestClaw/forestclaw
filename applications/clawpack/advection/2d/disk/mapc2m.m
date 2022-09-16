function [xp,yp,zp] = mapc2m(xc,yc)

parms = read_vars();

if parms.example == 0
    map = 'pillowdisk';
elseif parms.example == 1
    map = 'pillowdisk5';
end

alpha = parms.alpha;

shift = [1,1,0];
scale = [1,1,1];

switch map
    case 'pillowdisk'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        [xp,yp,~]= mapc2m_pillowdisk(xc,yc);
        xp = scale(1)*xp + shift(1);
        yp = scale(2)*yp + shift(2);
        
    case 'pillowdisk5'
        % (xc,yc) in [0,1]x[0,1]
        [xp,yp,~] = mapc2m_pillowdisk5(xc,yc,alpha);
        xp = scale(1)*xp+shift(1);
        yp = scale(2)*yp+shift(2);
end
zp = 0*xp;

end