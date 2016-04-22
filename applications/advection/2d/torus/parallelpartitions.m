function pp = parallelpartitions(q)
% PARALLELPARTITIONS returns a structure used to plot parallel partitions
%
% PP = PARALLELPARTITIONS(Q) returns basic information about the 

% pp.337;      % Good for four colors (pink, blue, green, ...)


pp.seed = 337;    % Random seed for processor colormap.
pp.npmax = 12;    % Number of processors

pp.qmin = 0.1;
pp.qmax = 0.9;
% pp.qmin = 2;
% pp.qmax = 4;

cm = [166	206	227;
    31	120	180;
    178	223	138;
    51	160	44;
    251	154	153;
    227	26	28;
    253	191	111;
    255	127	0;
    202	178	214;
    106	61	154;
    255	255	153;
    177	89	40;
    231, 41, 138;
    217, 95,2;
    27, 158,119;
    117, 112,179]/255;

pp.colormap_pp = cm(randperm(12),:);

pp.qcolors = q;
pp.colormap = colormap(yrbcolormap);  % Color map for q portion.

end