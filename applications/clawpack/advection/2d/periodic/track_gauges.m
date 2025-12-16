function [track_out,tcurr_out] = track_gauges(tcurr,gauges2track, format)
%
% track_gauges tracks the gauge locations for moving gauges over interval
% tspan=[t0,t1]
%
% track_gauges() will track all gauges using gauge locations from
% gaugeXXXXX.txt.  The full path will be shown, with dot for the current
% time. 
%
% track_gauges(gauges2track,tcurr,format) The format can be specified to 
% make sure that gauges are plotted at suitable z-level and that they will 
% appear  on top of all levels. Available formats are either 'geoclaw' or 
% 'forestclaw' (default).
%
% [track_hdl,tcurr_hdl] = track_gauges(...) returns a handle to the gauge 
% path and symbol at the current time. 
% 
% See also add_gauges, add_regions. 


if nargin < 3
    format = 'ForestClaw';
    if nargin < 2
        gauges2track = 'all';
    end
end

o = findobj(gcf,'Tag','gauge');
if (~isempty(o))
    delete(o);
end

o = findobj(gcf,'Tag','gauge_track');
if (~isempty(o))
    delete(o);
end


use_forestclaw = strcmpi(format,'forestclaw');

gauges = read_gauge_data();

np = get(gca,'NextPlot');
set(gca,'NextPlot','add');

% Set z-levels appropriately for each code. 
if (use_forestclaw)
    % ForestClaw plots coarser levels above finer levels so that 
    % patch boundaries show up nicely
    zmax = 0;
    zl = [-20,0];
else
    % AMRClaw plots finer levels on top of coarser levels.
    zmax = 20;
    zl = [0,20];
end

num_gauges = length(gauges);
gauge_ids = zeros(num_gauges,2);
for i = 1:num_gauges
    gid = gauges(i).id;
    gauge_ids(i,1) = gid;
    if strcmp(gauges2track,'all') == 1
        gauge_ids(i,2) = true;    
    else
        loc = find(gauges2track==gid);
        gauge_ids(i,2) = ~isempty(loc);
    end
end

outdir = './';
track_hdl = zeros(num_gauges,1);
tcurr_hdl = zeros(num_gauges,1);
for n = 1:num_gauges
    if ~gauge_ids(n,2)
        continue;
    end
    g = gauges(n);

    fprintf('Plotting gauges %d\n',g.id);
    [xvec,yvec,zvec,xcurr,ycurr,zcurr] = ...
        read_gauge_track(tcurr,g,outdir);

    % Plot location at the current time
    zp = zmax;

    % Plot gauge track
    hg1 = plot3(xvec,yvec,zp + 0*xvec,'k-','linewidth',2);
    set(hg1,'Tag','gauge_track')
    set(hg1,'userdata',g)

    % Plot location of gauge at current time
    hg2 = plot3(xcurr,ycurr,zp,'m.','linewidth',3,'markersize',50);
    set(gca,'zlim',zl);
    view(2);

    % set(gca,'zlimmode','auto');
    set(hg2,'Tag','gauge');
    g.h = text(xcurr,ycurr,zp,sprintf('%d',g.id),'fontsize',11,'color','k','fontweight','bold');
    set(g.h,'HorizontalAlignment','center');
    % set(h,'backgroundcolor','none');
    set(hg2,'userdata',g);
    tcurr_hdl(n) = hg2;   
    
end

zl = zlim;
if (use_forestclaw)
    set(gca,'zlim',[min(zl),zmax]);
else
    set(gca,'zlim',[0,zmax]);
end

% Restore hold status
set(gca,'NextPlot',np);
set(gca,'userdata',{tcurr_hdl,track_hdl});

if (nargout > 0)
    tcurr_out = tcurr_hdl;
    track_out = track_hdl;
end

end

function [xvec,yvec,zvec,xcurr,ycurr,zcurr] = ...
    read_gauge_track(tcurr,g,outdir)

gname = sprintf('%s/gauge%05d.txt',outdir,g.id);
if (exist(gname,'file'))

    % Skip '2' header lines
    tseries = importdata(gname,' ',2);

    % Index locations of (x,y,z) and tseries
    x_idx = 2;
    y_idx = 3;
    z_idx = 4;
    t_idx = 5;


    % Time series, and vector x,y,z values
    tvec = tseries.data(:,t_idx); 
    xvec = tseries.data(:,x_idx);
    yvec = tseries.data(:,y_idx);
    zvec = tseries.data(:,z_idx);
    
    L = length(tseries.data);
    
    if tcurr < tvec(1) || tcurr > tvec(end)
        xcurr = nan;
        ycurr = nan;
        zcurr = nan;
        return
    end

    % Find last entry in tvec that satifies 'tvec < tcurr'.
    tloc_min = find(tvec <= tcurr,1,'last');
    if tloc_min < L        
        % tvec(tloc_min) <= tcurr < tvec(tloc_min+1).
        tloc_max = tloc_min + 1;
        t0 = tvec(tloc_min);
        t1 = tvec(tloc_max);

        s = (tcurr - t0)/(t1-t0);

        x0 = tseries.data(tloc_min,x_idx);
        x1 = tseries.data(tloc_max,x_idx);
        xcurr = x0 + s*(x1-x0);

        y0 = tseries.data(tloc_min,y_idx);
        y1 = tseries.data(tloc_max,y_idx);
        ycurr = y0 + s*(y1-y0);

        z0 = tseries.data(tloc_min,z_idx);
        z1 = tseries.data(tloc_max,z_idx);
        zcurr = z0 + s*(z1-z0);
    
    else
        % tcurr is last entry  in tvec. 
        xcurr = tseries.data(tloc_min,x_idx);
        ycurr = tseries.data(tloc_min,y_idx);
        zcurr = tseries.data(tloc_min,z_idx);
    end
    %{    
    if ischar(plot_var)
        qval = feval(plot_var,tseries.data(:,t_idx+1:end));
    else
        qval = tseries.data(:,t_idx+plot_var);
    end
    %}

    x = tseries.data(:,x_idx);
    y = tseries.data(:,y_idx);
    z = tseries.data(:,z_idx);

else
    error('File %s does not exist\n',gname);
end

end



function gauges = read_gauge_data_XXX()

if (~exist('gauges.data','file'))
    fprintf('File gauges.data does not exist.  No gauges will be plotted.\n');
    gauges = [];
    return
end

fid = fopen('gauges.data','r');
for i = 1:5
    % Read first five lines of comments
    fgetl(fid);
end

gtype = struct('id',[],'x',[],'y',[],'t0',[],'t1',[]);

fgetl(fid);  % blank line
l = fgetl(fid);  % Dimension
dim = sscanf(l,'%d',1);
l = fgetl(fid);  % Get number of gauges
num_gauges = sscanf(l,'%d',1);
gauges(1:num_gauges) = gtype;
for n = 1:num_gauges
    l = fgetl(fid);
    data = sscanf(l,'%d %e %e %e %d',Inf);
    g = gtype;
    g.id = data(1);
    g.x = data(2);
    g.y = data(3);
    if dim == 2
        g.t0 = data(4);
        g.t1 = data(5);
    else
        g.z = data(4);
        g.t0 = data(5);
        g.t1 = data(6);
    end

    gauges(n) = g;
end

end
