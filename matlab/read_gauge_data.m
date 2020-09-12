

function gauges = read_gauge_data()
%
% read_gauge_data reads data in 'gauges.data' file and returns data in a
% struct.
%
% read_gauge_data() reads data in a file 'gauges.data'.  Data is assumed to
% be in the current directory. 
% 
% See also add_gauges, add_regions.

if (~exist('gauges.data','file'))
    fprintf('File gauges.data does not exist.  No gauges will be plotted.\n');
    return
end

gtype = struct('id',[],'longitude',[],'latitude',[],'t0',[],'t1',[]);

fid = fopen('gauges.data','r');
for i = 1:5
    % Read first five lines of comments
    fgetl(fid);
end

fgetl(fid);      % blank line
l = fgetl(fid);  % Get number of gauges
num_gauges = sscanf(l,'%d',1);
gauges(1:num_gauges) = gtype;
for n = 1:num_gauges    
    l = fgetl(fid);
    data = sscanf(l,'%d %e %e %e %d',Inf);
    g = gtype;
    g.id = data(1);
    g.longitude = data(2);
    g.latitude = data(3);
    g.t0 = data(4);
    g.t1 = data(5);
    gauges(n) = g;
end

end