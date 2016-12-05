function eta = bathy(data)
% read q data:
h = data(:,1);
eta = data(:,4);
eta(h <= 1e-3) = nan;

end
