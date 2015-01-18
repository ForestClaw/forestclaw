function zm = mountain_height(xp)

mt = load('mountain.dat');

zm = pchip(mt(:,1),mt(:,2),xp);

end

