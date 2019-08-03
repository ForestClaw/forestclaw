function h = plotframe2ez(amrdata,mq,plotstyle,map1d)

%
% PLOTFRAME2EZ plots data from a single Clawpack output file.
%
%    PLOTFRAME2EZ 

%    There one basic type of plot available.  These are
%
%         Surface plots : 2d plots of the data.  This data may be viewed as
%         a manifold, if desired.
%
%         Scatter plots : Spherically symmetric data is plotted as a
%         function of a 1d variable, usually the distance from a fixed point
%         (x0,y0).  The result is a line-type plot of points representing the
%         data.
%
%    See SETPLOT, PLOTCLAW2, MANIFOLD.

if isempty(amrdata)
    fprintf('\n');
    fprintf('Frame %d (%s) does not exist\n',Frame,outputflag);
    fprintf('\n');
    return;
end


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
