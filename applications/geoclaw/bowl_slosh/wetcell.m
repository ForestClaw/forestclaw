function eta = wetcell(data)
% read q data:
mbathy = 4;
mheight = 1;
h = data(:,mheight);
b = data(:,mbathy);
eta = zeros(size(h));

eta = data(:,4);

dry_tol = 0.001;
eta(h <= dry_tol) = nan;
end