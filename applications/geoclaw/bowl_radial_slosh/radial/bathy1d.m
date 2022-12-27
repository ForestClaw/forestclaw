function eta = bathy1d(data)
% Needs to be modified for 1d problem.
mbathy = 4;
mheight = 1;
h = data(:,mheight);
b = zeros(size(h));   % Needs to be modifed
eta = zeros(size(h));
m = b > 0;
eta(b>0) = nan;

eta(b<=0) = h(b<=0) + b(b<=0);
end