function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 


pp.seed = 310;    % okay
pp.seed = 312;   % Random seed for processor colormap.
pp.npmax = 5;    % Number of processors

tol = 1e-2;
pp.qmin = 0.01;
pp.qmax = 1.1;
pp.plotq = true;

pp.qcolors = q;
pp.colormap = colormap(yrbcolormap);  % Color map for q portion.

end