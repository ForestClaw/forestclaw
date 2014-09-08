function [xp,yp,zp] = mapc2m_squareddisk(xc,yc)

blockno = getblocknumber();

alpha = 0.5;

if (blockno == 2)
    xp = (2*xc - 1)*alpha/sqrt(2);
    yp = (2*yc - 1)*alpha/sqrt(2);
else
    switch blockno
        case 0
            xc1 = xc;
            yc1 = 1-yc;
            [xp,yp] = squareddisk_help(alpha,xc1,yc1);
            yp = -yp;
        case 1
             xc1 = yc;
             yc1 = 1-xc;
             [yp,xp] = squareddisk_help(alpha,xc1,yc1);
             xp = -xp;            
        case 3
             xc1 = yc;
             yc1 = xc;
             [yp,xp] = squareddisk_help(alpha,xc1,yc1);
        case 4
             xc1 = xc;
             yc1 = yc;
             [xp,yp] = squareddisk_help(alpha,xc1,yc1);
        otherwise
            str = sprintf('blockno = %d\n',blockno);
            error(str);
    end
end

zp = 0*xp;
    
end

function [x,y] = squareddisk_help(alpha,xi,eta)

% Original mapping
R = alpha.^(1-eta);

% This appears to lead to more uniform cells...
% R = (1-alpha)*eta + alpha;

tan_xi = tan(0.5*pi*(xi - 0.5));
xi_prime = 2*(1 - eta).*(xi-0.5)+eta.*tan_xi;
y = R./sqrt(1 + eta.*tan_xi.^2 + (1-eta));
x = y.*xi_prime;
end




