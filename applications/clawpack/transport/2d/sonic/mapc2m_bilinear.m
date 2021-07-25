function [xp,yp,zp] = mapc2m_bilinear(xc,yc,center)

blockno = getblocknumber();

quad(1,1,1:2) = [0;0];
quad(2,1,1:2) = [1;0];
quad(1,2,1:2) = [0;1];
quad(2,2,1:2) = [1;1];

switch blockno
    case 0
        s = [-1;-1];        
    case 1
        s = [0; -1];
    case 2
        s = [-1; 0];
    case 3
        s = [0;0];
end

for i = 1:2
    for j = 1:2
        quad(i,j,1) = quad(i,j,1) + s(1);
        quad(i,j,2) = quad(i,j,2) + s(2);
    end
end

switch blockno
    case 0
        quad(2,2,:) = center(:);
    case 1
        quad(1,2,:) = center(:);
    case 2
        quad(2,1,:) = center(:);
    case 3
        quad(1,1,:) = center(:);
end

[m,n] = size(xc);

xp = zeros(size(xc));
yp = zeros(size(yc));
for i = 1:m
    for j = 1:n
        for k = 1:2
            a00(k) = quad(1,1,k);
            a01(k) = quad(2,1,k) - quad(1,1,k);
            a10(k) = quad(1,2,k) - quad(1,1,k);
            a11(k) = quad(2,2,k) - quad(2,1,k) - quad(1,2,k) + quad(1,1,k);
        end
        xp(i,j) = a00(1) + a01(1)*xc(i,j) + a10(1)*yc(i,j) + a11(1)*xc(i,j)*yc(i,j);
        yp(i,j) = a00(2) + a01(2)*xc(i,j) + a10(2)*yc(i,j) + a11(2)*xc(i,j)*yc(i,j);
    end    
end
zp = zeros(size(xc));

end
            
            