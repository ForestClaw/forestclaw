function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 


pp.seed = 337;   % Random seed for processor colormap.
pp.npmax = 5;    % Number of processors

pp.qmin = 0.18;
pp.qmax = 3.4169;
pp.qmin = 1e-16;
pp.qmax = 3e-4;
pp.plotq = true;

pp.qcolors = q;
pp.colormap = colormap(yrbcolormap);  % Color map for q portion.

end