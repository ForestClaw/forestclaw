function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 
% patches should be colored.  Q is the data to be used for color range.
%

pp.seed = 337;   % Random seed for processor colormap.
pp.npmax = 4;    % Number of processors
% pp.pp_colormap = [1 0.6 0.6; 0.6 0.6 1];

pp.qmin = 0.205;
pp.qmax = 0.25;

pp.colormap = colormap(jet(64));  % Color map for q portion.

pp.plotq = true;

pp.qcolors = [];

if (nargin == 0)
    return
end

qcolors = q;

% Map q linearly into colormap; constrain entire color map into
% range [qmin, qmax].  Values outside of this range will use
% the processor color.

pp.qcolors = qcolors;


end