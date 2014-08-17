function [xp,yp,zp] = mapc2m_squareddisk(xc,yc)

%       subroutine mapc2m_squareddisk(xc,yc,xp,yp,zp)
%       implicit none
% 
%       double precision xc,yc,xp,yp,zp
%       double precision user_double(0:15)
%       double precision half_length
%       double precision R2sqrbyR1, R1byR2
%       double precision xc1, yc1
%       integer blockno
% 
%       blockno = get_block()
% 
%       if (blockno .eq. 2) then
%          half_length = user_double(2);
%          xp = (2*xc - 1)*half_length;
%          yp = (2*yc - 1)*half_length;
%       else
%          R2sqrbyR1 = double_user(0)
%          R1byR2 = user_double(1)
% 
%          if (blockno .eq. 0) then
%             xc1 = xc
%             yc1 = 1.d0-yc
%             call squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1,xp,yp)
%             yp = -yp
%          elseif (blockno .eq. 1) then
%             xc1 = yc
%             yc1 = 1.d0-xc
%             call squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1,xp,yp)
%             xp = -xp
%          elseif (blockno .eq. 3) then
%             xc1 = yc
%             yc1 = xc
%             call squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1,xp,yp)
%          elseif (blockno .eq. 4) then
%             xc1 = xc
%             yc1 = yc
%             call squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1,xp,yp)
%          endif
%       end
% 
% 
%       subroutine squareddisk_help(R2sqrbyR1, R1byR2,
%      &      xi, eta, x, y)
%       implicit none
% 
%       double precision R2sqrbyR1, R1byR2
%       double precision R, tan_xi, xi_prime
% 
%       double precision pi
% 
%       pi = 3.1415926535897932384626433d0
% 
%       R = R2sqrbyR1*R1byR2**(1.d0 + eta)
%       tan_xi = tan(0.5d0*pi*(xi - 0.5d0))
%       xi_prime = (1.d0 - eta)*2.d0*(xi-0.5d0)+eta*tan_xi
% 
%       y = R/sqrt(1.d0 + eta*tan_xi**2 + (1.d0-eta))
%       x = y*xi_prime
% 
% 
%       end

blockno = getblocknumber();

R1 = 1;
R2 = 0.5;

R2sqrbyR1 = R2^2/R1;
R1byR2 = R1/R2;
half_length = R2/sqrt(2);

if (blockno == 2)
    xp = (2*xc - 1)*half_length;
    yp = (2*yc - 1)*half_length;
else
    switch blockno
        case 0
            xc1 = xc;
            yc1 = 1-yc;
            [xp,yp] = squareddisk_help(R2sqrbyR1,R1byR2,xc1,yc1);
            yp = -yp;
        case 1
             xc1 = yc;
             yc1 = 1-xc;
             [yp,xp] = squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1);
             xp = -xp;            
        case 3
             xc1 = yc;
             yc1 = xc;
             [yp,xp] = squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1);
        case 4
             xc1 = xc;
             yc1 = yc;
             [xp,yp] = squareddisk_help(R2sqrbyR1, R1byR2,xc1,yc1);
        otherwise
            str = sprintf('blockno = %d\n',blockno);
            error(str);
    end
end

zp = 0*xp;
    
end

function [x,y] = squareddisk_help(R2sqrbyR1, R1byR2,xi,eta)

R = R2sqrbyR1*R1byR2.^(1 + eta);
tan_xi = tan(0.5*pi*(xi - 0.5));
xi_prime = 2*(1 - eta).*(xi-0.5)+eta.*tan_xi;
y = R./sqrt(1 + eta.*tan_xi.^2 + (1-eta));
x = y.*xi_prime;
end




