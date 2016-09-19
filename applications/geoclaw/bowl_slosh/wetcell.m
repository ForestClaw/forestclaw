function eta = wetcell(data)
% read q data:
mbathy = 4;
mheight = 1;
h = data(:,mheight);
b = data(:,mbathy);
eta = zeros(size(h));

dry_tol = 0.001;
eta(h>dry_tol) = h(h>dry_tol)+b(h>dry_tol);
eta(h<=dry_tol) = nan;
end