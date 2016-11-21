function eta = bathy(data)
% read q data:
h = data(:,1);
eta = data(:,4);
% eta = zeros(size(h));
eta(h <= 0) = nan;

end