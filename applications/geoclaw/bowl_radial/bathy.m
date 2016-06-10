function plotobject = bathy(data)
% read q data:
mbathy = 4;
mheight = 1;
plotobject = zeros(size(data(:,1)));
plotobject(data(:,mbathy)>0) = -1;
plotobject(data(:,mbathy)<=0) = data(data(:,mbathy)<=0,mheight) ...
    + data(data(:,mbathy)<=0,mbathy);
end