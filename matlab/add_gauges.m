function ghout = add_gauges()

o = findobj(gcf,'Tag','gauge');
if (~isempty(o))
    delete(o);
end

np = get(gca,'NextPlot');

if (~exist('gauges.data'))
    fprintf('File gauges.data does not exist.  No gauges will be plotted.\n');
    return
end

fid = fopen('gauges.data','r');
for i = 1:5,
    % Read first five lines of comments
    l = fgetl(fid);
end
l = fgetl(fid);  % blank line
l = fgetl(fid);  % Get number of gauges
num_gauges = sscanf(l,'%d',1);
gauge_handles = zeros(num_gauges,1);
for n = 1:num_gauges,
    l = fgetl(fid);
    data = sscanf(l,'%d %e %e %e %d',Inf);
    hg = plot(data(2),data(3),'x','linewidth',3,'markersize',8);
    set(hg,'Tag','gauge');
    hold on;
    h = text(data(2),data(3),sprintf('  %d',data(1)),'fontsize',16);
    set(h,'HorizontalAlignment','right');
end

set(gca,'NextPlot',np);

if (nargout > 0)
    ghout = gauge_handles;
end

end