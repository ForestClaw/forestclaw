function eta = wetcell(data)
% read q data:
mheight = 1;
meta = 4;

h = data(:,mheight);

eta = data(:,meta);

dry_tol = 0.001;
eta(h <= dry_tol) = nan;

end