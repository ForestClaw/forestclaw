function plot_gauges(gauges2plot,plot_var, outdir, show_levels)
% plot_gauges - Create plots of data stored in `gaugeXXXXX.txt` files.
% 
% 
% Syntax
%   plot_gauges()
%   plot_gauges(gauges2plot)
%   plot_gauges(gauges2plot,plot_var)
%   plot_gauges(gauges2plot,plot_var,outdir)
%   plot_gauges(gauges2plot,plot_var,outdir,show_levels)
% 
% Input Arguments
%   gauges2plot - Array of gauges to plot.  Supply a scalar, 
%   an array of gauge IDs, or string `all` to plot all gauges. Default is
%   `all`. 
% 
% plot_gauges(gauges2plot) plots variable 'plot_var', where `plot_var`
% corresponds to either an index in q-state, or a string referring to a
% function to call to compute derived quantity. 
% 
% plot_gauges(...,outdir) specifies a directory where gauge files are
% stored.  
% 
% Examples : 
%   To plot "energy" for all gauges
% 
%        plot_gauges('all');
%
%   or just 
% 
%        plot_gauges(4);
% 
%   To plot pressure at all gauges
% 
%        plot_gauges('pressure');
% 
%   To plot pressure at a subset o gauges
% 
%        plot_gauges('pressure',[111,232,353]);
% 

if nargin < 4
    show_levels = false;
    if nargin < 3
        plot_var = 1;
        if nargin < 2
            outdir = './';
            if nargin < 1
                gauges2plot = 'all';
            end
        end
    end
end

close all

% ----------------------------------
% Read 'gauges.data' file
% ----------------------------------
gauges = read_gauge_data();

%{
if ischar(gnum)
    gauge_list = [gdata.id];
elseif isscalar(gnum)
    gauge_list = gnum;
end    
%}

num_gauges = length(gauges);
gauge_ids = zeros(num_gauges,2);
for i = 1:num_gauges
    gid = gauges(i).id;
    gauge_ids(i,1) = gid;
    if strcmp(gauges2plot,'all') == 1
        gauge_ids(i,2) = true;    
    else
        loc = find(gauges2plot==gid);
        gauge_ids(i,2) = ~isempty(loc);
    end
end


if ischar(plot_var)
    ystr = plot_var;
elseif isnumeric(plot_var)
    var_idx = plot_var + 5;   % Could be an aux variable.
%{
    if plot_var <= meqn
        ystr = 'q';
        var_idx = plot_var;
    elseif meqn < plot_var && plot_var <= meqn + maux
        ystr = 'aux';
        var_idx = plot_var-meqn;
    end
%}
else
    error('Invalid variable specifed.  Specify a string or variable in [1-%d]',meqn+maux);
end

clist = {'r','b','g','b','c','m'};

for i = 1:num_gauges    
    if ~gauge_ids(i,2)
        continue;
    end
    g = gauges(i);

    figure(100 + i);
    clf;
    hold on;
    
    gname = sprintf('%s/gauge%05d.txt',outdir,g.id);
    if (exist(gname,'file'))

        % Skip 2 lines in the header
        tseries = importdata(gname,' ',2);

        t_idx = 5;
        if ~isfield(tseries,'data')
            fprintf('Gauge file for gauge ID %d does not have any data.\n',g.id);
        else
            tvec = tseries.data(:,t_idx);

            if ischar(plot_var)
                qval = feval(plot_var,tseries.data(:,t_idx+1:end));
            else
                qval = tseries.data(:,t_idx+plot_var);
            end
            plot(tvec,qval,'k.-','linewidth',1,'markersize',8);
            if show_levels
                c = tseries.data(:,1);  % levels
                lvals = sort(unique(c));
                hdl = [];
                lstr = {};
                for i = 1:length(lvals)
                    l = lvals(i);
                    m = c == l;
                    hdl(i) = plot(t(m),qval(m),[clist{i},'.'],'markersize',10);
                    lstr{i} = sprintf("Level %d",l);
                end
            end
        end
    else
        error('File %s does not exist\n',gname);
    end
                              
    title(sprintf('Gauge %d',g.id),'fontsize',18);
    xlabel('t (seconds)','fontsize',16);
    % ylabel(ystr,'fontsize',16);
    set(gca,'fontsize',16);

    if show_levels
        legend(hdl,lstr);
    end
           
    hold off;
    shg
       
end

end

function gauges = read_gauge_data()
%
% read_gauge_data() reads data in a file 'gauges.data'.  Data is assumed to
% be in the current directory. 
% 
% See also add_gauges, add_regions.

if (~exist('gauges.data','file'))
    fprintf('File gauges.data does not exist.  No gauges will be plotted.\n');
    return
end

gtype = struct('id',[],'dim',[],'x',[],'y',[],'t0',[],'t1',[]);

fid = fopen('gauges.data','r');
for i = 1:5
    % Read first five lines of comments
    fgetl(fid);
end

fgetl(fid);      % blank line
l = fgetl(fid);  % Get number of gauges
dim = sscanf(l,'%d',1);
l = fgetl(fid);  % Get number of gauges
num_gauges = sscanf(l,'%d',1);
gauges(1:num_gauges) = gtype;
for n = 1:num_gauges    
    l = fgetl(fid);
    data = sscanf(l,'%d %e %e %e %d',Inf);
    g = gtype;
    g.dim = dim;
    g.id = data(1);
    g.x = data(2);
    g.y = data(3);
    if dim == 2
        g.t0 = data(4);
        g.t1 = data(5);
    elseif dim == 3
        g.z = data(4);
        g.t0 = data(5);
        g.t1 = data(6);
    end
    gauges(n) = g;
end

end
