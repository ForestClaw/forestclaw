function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the



pp.seed = 337;   % Random seed for processor colormap.
pp.npmax = 8;    % Number of processors

pp.qmin = 0.01;
pp.qmax = 1.01;
pp.qmin = 2;
pp.qmax = 4;

pp.qcolors = q;
pp.colormap = colormap(yrbcolormap);  % Color map for q portion.

end
