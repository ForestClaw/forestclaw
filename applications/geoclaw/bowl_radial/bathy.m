function eta = bathy(data)
% read q data:
mbathy = 4;
mheight = 1;
h = data(:,1);
b = data(:,4);
eta = zeros(size(h));
m = b > 0;
eta(b>0) = nan;

eta(b<=0) = h(b<=0) + b(b<=0);
end