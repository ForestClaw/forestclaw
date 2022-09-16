function [xp,yp,zp] = mapc2m_fivepatch(xc,yc,alpha)
blockno = getblocknumber();

[m,n] = size(xc);

if (nargin < 3)
    alpha = 0.5;
end

if (blockno == 2)
    xp = (2*xc - 1)*alpha;
    yp = (2*yc - 1)*alpha;
else
    switch blockno
        case 0
            xc1 = xc;
            yc1 = 1-yc;
            [xp,yp] = bilinear_help(alpha,xc1(:),yc1(:));
            yp = -yp;
        case 1
             xc1 = yc;
             yc1 = 1-xc;
             [yp,xp] = bilinear_help(alpha,xc1(:),yc1(:));                          
             xp = -xp;            
        case 3
             xc1 = yc;
             yc1 = xc;
             [yp,xp] = bilinear_help(alpha,xc1(:),yc1(:));  
        case 4
             xc1 = xc;
             yc1 = yc;
             [xp,yp] = bilinear_help(alpha,xc1(:),yc1(:));
        otherwise
            str = sprintf('blockno = %d\n',blockno);
            error(str);
    end
end
xp = reshape(xp,m,n);
yp = reshape(yp,m,n);
zp = 0*xp;
    
end

function [xb,yb] = bilinear_help(alpha,xi,eta)

% xpc = [-alpha, alpha 1 -1];
% ypc = [alpha alpha 1 1];

xpc = [-alpha, -1, 1, alpha];
ypc = [alpha 1 1 alpha];



a = [xpc(1); ypc(1)];
u1 = [xpc(4) - xpc(1); ypc(4) - ypc(1)];
v1 = [xpc(2) - xpc(1); ypc(2) - ypc(1)];
v2 = [xpc(3) - xpc(4); ypc(3) - ypc(4)];

len = length(xi(:));

xb = zeros(size(xi));
yb = zeros(size(eta));
% for k = 1:len,
%   pt = a + u1*xi(k) + v1*eta(k) + (v2-v1)*xi(k)*eta(k);
%   xb(k) = pt(1);
%   yb(k) = pt(2);
% end
xb = a(1) + u1(1)*xi + v1(1)*eta + (v2(1)-v1(1))*xi.*eta;
yb = a(2) + u1(2)*xi + v1(2)*eta + (v2(2)-v1(2))*xi.*eta;

end




