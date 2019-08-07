function [xp,yp,zp] = mapc2m(xc,yc)

[example,mapping,ic,alpha,center] = read_vars();

map_list = {'identity','cart','fivepatch','bilinear'};
map = map_list{mapping+1};

shift = [0.5,0.5,0];
scale = 2;

switch map
    case 'identity'
        % Uses [ax,bx] x [ay,by] set in fclaw2d_defaults.ini
        % This is what is stored in the fort.q files.
        % (xc,yc) in [0,1]x[0,1]
        s = 0.005;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
    
        xp = xc1;
        yp = yc1;
     
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0.005;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,~] = mapc2m_cart(xc1,yc1);
                        
        xp = xp/scale + shift(1);
        yp = yp/scale + shift(2);
        
    case 'fivepatch'
        alpha = 0.4;
        [xp,yp,~] = mapc2m_fivepatch(xc,yc,alpha);
        b = getblocknumber();
        s = 0.005;
        switch b
            case 0
                yp = yp - s;
            case 1
                xp = xp - s;
            case 3
                xp = xp + s;
            case 4
                yp = yp + s;
        end
        xp = xp/scale + shift(1);
        yp = yp/scale + shift(2);
        
    case 'bilinear'
        center = load('center.dat');
        [xp,yp,~] = mapc2m_bilinear(xc,yc,center);
        xp = xp/scale + shift(1);
        yp = yp/scale + shift(2);
end
zp = 0*xp;




end