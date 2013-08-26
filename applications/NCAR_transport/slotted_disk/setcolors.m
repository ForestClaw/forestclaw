function setcolors(p,x,y,z,q)
%
% This colormap sets shows partition colors.  Use in
% conjunction with MULTICOLORMAP.
% ToDo : Add this to plotframe to have a general routine
% for specifying partition colors, e.g.
%
%     ShowParallelPartitions = 1
%     PartitionField = 2     % mq value for partitions
%

% Use rank number to map directly into colormap.
plot_type = 'mesh';

switch plot_type
case 'partition'
cdata_idx = round(q) + 1;    % Integer partition values
set(p,'CData',cdata_idx);

set(p,'CDataMapping','direct');  % Scale into color map (set in afterframe).
set(p,'FaceColor','flat');       % Single color per cell

case 'mesh'
cdata_idx = 0*q + 1;
set(p,'CData',cdata_idx);
set(p,'CDataMapping','direct');  % Scale into color map (set in afterframe).
set(p,'FaceColor','flat');       % Single color per cell

end
