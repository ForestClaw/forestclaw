OutputFlag = 'ForestClaw';         % default value.
OutputDir = './';            % Default (reassign them here anyway)

ForestClaw = 1;     % Plot using ForestClaw preferences.

mq = 1;
UserVariable = 1;
UserVariableFile = 'surface_eta';
MaxFrames = 1000;
MaxLevels = 30;
Manifold = 0;
ReadBlockNumber = 1;
PlotData =  ones(1,MaxLevels);
PlotGrid =  zeros(1,MaxLevels);
PlotGridEdges = zeros(1,MaxLevels);

PlotType = 1;  % Scatter plot

MappedGrid = 0;

% for contour plots (PlotType==2):
ContourValues = [];

% for scatter plot (PlotType==4):
x0 = 0;
y0 = 0;
UserMap1d = 0;
ScatterStyle = setplotstyle('b-','r*','cs','gp','bo','k*');

PlotParallelPartitions=0;
