function h = plotframe2ez(amrdata,mq,plotstyle,map1d)
% 
% plotframe2ez plots 2d data from a specified data file.
%
%    plotframe2ez(AMRDATA,mq,plotstyle,map1d) takes a 1d slice, or 
%    some other mapping R2 --> R1 of data in each AMR grid and 
%    adds it to the current plot.  
%
%    The main purpose of this is routine is to compare slices from 
%    two different runs.  
%
%    H = PLOTFRAME2EZ(AMRDATA,mq,plotstyle,map1d) returns a vector of
%    handles to the line or scatter plots created from each grid.
%
%    NOTE : This should only work with PlotType == 4.  A function 'MAP1D'
%    should be written that maps data in R^2 to data in R. 
%
%    Example : 
%
%    % This example compare 2d data stored in 'fort_compare' to data
%    % plotted in the current call to plotclaw2.
%    %
%    % The following should go in afterframe.m
%    
%    [amrdata_ref,tref] = readamrdata_forestclaw(2,Frame,'./fort_compare/');  
%    lstyle = {'ro-','go-','bo-','mo-'};
%    href = plotframe2ez(amrdata_ref,mq,lstyle,@map1d);
%
%    
%
%    See SETPLOT, PLOTFRAME1EZ.

% -------------
% Scatter plots
% -------------

pstyle = plotstyle;
maxlevels = 30;

[linestyle,linecolors,markerstyle] = get_plotstyle(pstyle,maxlevels);

qmin = [];
qmax = [];
ncells = zeros(maxlevels,1);

%=============================================
% MAIN LOOP ON GRIDS FOR THIS FRAME:
%=============================================

forestclaw = true;
UserVariable = 0;

ngrids = length(amrdata);  % length of structure array
h = zeros(ngrids,1);
for ng = 1:ngrids
    
    % gridno = amrdata(ng).gridno;
    % blockno = amrdata(ng).blockno;   % == 0 if there is only one block
    level = amrdata(ng).level;
    if (forestclaw)
        level = level + 1;   % ForestClaw levels start at 0
    end
    
    mx = amrdata(ng).mx;
    my = amrdata(ng).my;
    
    xlow = amrdata(ng).xlow;
    ylow = amrdata(ng).ylow;
    
    dx = amrdata(ng).dx;
    dy = amrdata(ng).dy;
    
    xedge = xlow + (0:mx)*dx;
    yedge = ylow + (0:my)*dy;
    
    xcenter = xedge(1:mx) + dx/2;
    ycenter = yedge(1:my) + dy/2;
    
    % for compatibility with old matlab41/plotframe2 convention:
    x = xcenter;
    y = ycenter;
    
    % read q data:
    data = amrdata(ng).data;
        
    data = data';
    
    if (UserVariable == 1)
        % User has supplied a function to convert original q variables to
        % the variable which is to be plotted, e.g. Mach number, entropy.
        qdata = feval(UserVariableFile,data);
        q = reshape(qdata,mx,my);
    else
        q = reshape(data(:,mq),mx,my);
    end
    
    amrdata(ng).q = q;
    
    qmesh = q';
    
    % minimum over all grids at this time, but not necessarily on slice
    % shown.
    qmin = min([qmin,min(q(:))]);
    qmax = max([qmax,max(q(:))]);
    
    % keep count of how many cells at this refinement level:
    if length(ncells) < level
        ncells(level) = 0;
    end
    ncells(level) = ncells(level) + mx*my;
    
    
    % 1d Line plots
    
    [xcm,ycm] = meshgrid(xcenter,ycenter);
    
    % Users should call mapc2p from inside of map1d.
    [rvec,qvec] = map1d(xcm,ycm,qmesh);
    cs = size(rvec,2);
    cq = size(qvec,2);
    if (cs > 1 || cq > 1)
        error(['plotframe2ez : map1d should return single column vectors ',...
            'for r or q']);
    end
 
    h(ng) = line(rvec,qvec,'Marker',markerstyle{level},'Color',...
        linecolors{level},'LineStyle',linestyle{level});
    
%    add_line2plot(rvec,qvec,level,markerstyle{level},...
%        linecolors{level},linestyle{level});
            
end % loop on ng (plot commands for each grid)
%=============================================

end
