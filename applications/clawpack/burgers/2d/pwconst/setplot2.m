ForestClaw = 1;
OutputFlag = 'forestclaw';
OutputDir = './';            % Default (reassign them here anyway)

PlotType = 1;                % type of plot to produce:
			     % 1 = pseudo-color (pcolor)
			     % 2 = contour
			     % 3 = Schlieren
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
Manifold = 0;
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 30;               % increase if using amrclaw with more levels
PlotData =  ones(MaxLevels,1);   % Data on refinement level k is plotted only if
			     % k'th component is nonzero
PlotGrid =  zeros(MaxLevels,1);   % Plot grid lines on each level?

PlotGridEdges =  ones(MaxLevels,1);  % Plot edges of patches of each grid at
                                 % this level?

ContourValues = [];

