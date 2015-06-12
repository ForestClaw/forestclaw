function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 



pp.seed = 337;   % Random seed for processor colormap.
pp.npmax = 4;    % Number of processors

pp.qcolors = [];

qcolors = q;
m = isnan(q);
qcolors(m) = 0.5;
qcolors(~m) = nan;  % Other values should be colored by processor color

% Map q linearly into colormap; constrain entire color map into
% range [qmin, qmax].  Set values that should use the partition
% colormap to nan;

pp.qmin = 0;
pp.qmax = 1;
pp.qcolors = qcolors;
pp.colormap = colormap(white(64));  % Color map for q portion.

end