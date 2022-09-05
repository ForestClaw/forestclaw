%  setplot3.m
%  called in plotclaw3 before plotting to set various parameters

OutputDir = './';
setviews  % set viewpoints so that view(xSlice), for example, can be used.
setopengl;

OutputFlag = 'forestclaw';

ForestClaw = 1;

PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r
			     % 5 = isosurface plot (at given levels)

mq = 1;                      % which component of q to plot
UserVariable = 1;            % set to 1 to specify a user-defined variable
UserVariableFile = 'pressure';      % name of m-file mapping data to q
MappedGrid = 1;              % set to 1 if mapc2p.m exists for nonuniform grid
Manifold = 0;
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 30;               % max number of AMR levels
ReadBlockNumber = 1;

PlotData =  ones(1,MaxLevels);       % Data on refinement level k is plotted only
			         % if k'th component is nonzero
PlotGrid =  ones(1,MaxLevels);       % Plot grid lines on each level?
PlotGridEdges =  zeros(1,MaxLevels);  % Plot edges of patches of each grid at
                                 % this level on slices?
PlotCubeEdges = zeros(1,MaxLevels);   % Plot edges of cube of refinement patch at
                                 % this level?


% ContourValues is a vector of contour lines that can be used with
% PlotType = 1,2,3.   Empty ==> no contour lines drawn

ContourValues = [];   % draw contour lines from 'afterframe.m'

% The next three parameters are vectors of x,y,z coordinates of 2d slices
% to be displayed for PlotType = 1,2,3.

% Note : Supply computational coordinates.  This is a bit of a disaster
% for multiblock domains ...
xSliceCoords = [];
ySliceCoords = [];
zSliceCoords = [];
parms = read_vars();
switch parms.mapping
    case 0
        zSliceCoords = [-0.5, 0, 0.5];
    case 1
        zSliceCoords = [0.25 0.5 0.75];   % [0.2 0.4 0.6 0.8];
    otherwise
        error("zSlicecoords not specified");
end

IsosurfValues    =  [1.05];        % Plot surfaces at q = surfValue(i).

IsosurfAlphas    =  [0.75];     % Transparency of each surface
                                          % (0=clear; 1=opaque)
                                          % NOTE: Your system must be able to
                                          % use the OpenGL Renderer.

IsosurfColors    = 'w';      % Colors for each surface.

PlotBlockColors = 0;
PlotParallelPartitions = 0;
