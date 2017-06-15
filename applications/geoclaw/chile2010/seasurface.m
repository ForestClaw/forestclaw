function eta = seasurface(data)
% read q data:
h = data(:,1);
eta = data(:,4);
b = eta - h;
bm = b < 0;
hm = h < 0;
eta(b > 0) = nan;

end
