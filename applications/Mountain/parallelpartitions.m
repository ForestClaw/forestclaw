function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 



pp.seed = 237;   % Random seed for processor colormap.
pp.npmax = 8;    % Number of processors

pp.qcolors = [];


% Map q linearly into colormap; constrain entire color map into
% range [qmin, qmax].  Colors outside of this range will use 
% the parallel partition colors. 

pp.qmin = -100;
pp.qmax = -50;
pp.colormap = colormap(jet(64));  % Color map for q portion.

end