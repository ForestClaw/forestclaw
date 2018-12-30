OutputFlag = 'ForestClaw';         % default value.
OutputDir = './';            % Default (reassign them here anyway)

ForestClaw = 1;     % Plot using ForestClaw preferences.

PlotType = 1;                % type of plot to produce:
			     % 1 = pseudo-color (pcolor)
			     % 2 = contour
			     % 3 = Schlieren
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
Manifold = 1;
ReadBlockNumber = 1;
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 30;
PlotData =  ones(1,MaxLevels);   % Data on refinement level k is plotted only if
			     % k'th component is nonzero
PlotGrid =  zeros(1,MaxLevels);   % Plot grid lines on each level?

PlotGridEdges = ones(1,MaxLevels);  % Plot edges of patches of each grid at
                                 % this level?

%---------------------------------

% Set to either a scalar, for automatic contours or a vector of contour levels.
ContourValues = [];

%---------------------------------

ShowUnderOverShoots = 1;
PlotParallelPartitions=0;
