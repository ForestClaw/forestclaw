function rhout = add_regions(t)

o = findobj(gcf,'Tag','region');
if (~isempty(o))
    delete(o);
end

np = get(gca,'NextPlot');

if (~exist('regions.data'))
    fprintf('File regions.data does not exist.  No regions will be plotted.\n');
    return
end

fid = fopen('regions.data','r');
for i = 1:5,
    % Read first five lines of comments
    l = fgetl(fid);
end
l = fgetl(fid);  % blank line
l = fgetl(fid);  % Get number of gauges
num_regions = sscanf(l,'%d',1);

c = {'r','b','y','c','w','g','k'};   % Colors for each region

region_handles = zeros(num_regions,1);
for n = 1:num_regions,
    l = fgetl(fid);
    data = sscanf(l,'%d %d %e %e %e %e %e %e',Inf);
    minlevel = data(1);
    maxlevel = data(2);
    t0 = data(3);
    t1 = data(4);
    x0 = data(5);
    x1 = data(6);
    y0 = data(7);
    y1 = data(8);
    xp = [x0 x1 x1 x0 x0];
    yp = [y0 y0 y1 y1 y0];
    if (t0 <= t && t <= t1)
        hg = plot(xp,yp,'linewidth',2,'color',c{mod(n-1,7)+1});
        set(hg,'Tag','region');
    end
    zl = zlim;
    set(gca,'zlim',[min(zl),0]);
    hold on;
end

set(gca,'NextPlot',np);

if (nargout > 0)
    rhout = region_handles;
end

end
