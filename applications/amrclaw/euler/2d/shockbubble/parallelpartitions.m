function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 



pp.seed = 337;   % Random seed for processor colormap.
pp.npmax = 5;    % Number of processors

pp.qmin = 0.18;
pp.qmax = 3.4169;
pp.plotq = false;

pp.qcolors = q;
pp.colormap = colormap(yrbcolormap);  % Color map for q portion.

end