function plot_gauges(plot_var)

if nargin < 1
    plot_var = 'v';
end

close all

plot_geo = true;
plot_forest = false;
plot_obs_data = true;

% ------------------------------------------
geo_outdir = './gauges_level_6_geo';
geo_outdir = './_from_fusion';

fc_outdir = './gauges_level_14/subcycle';
% outputdir = './';


% Scaling
tscale = 1/3600;     % Convert time to hours
if (plot_obs_data)
    tshift = 10*60;
end
vel_scale = 100;       % Convert velocity to cm/sec


% ----------------------------------
gdata = read_gauge_data();

t_idx = 2;
h_idx = 3;
hu_idx = 4;
hv_idx = 5;
eta_idx = 6;

switch plot_var
    case 'eta'
        pvidx = 1;
        ystr = 'Surface height';
    case 'u'
        pvidx = 2;
        ystr = 'u-velocity';
    case 'v'
        pvidx = 3;
        ystr = 'v-velocity';
    case 'speed'
        pvidx = 4;
        ystr = 'Speed';
end

ph = [];
lstr = {};
num_gauges = length(gdata);
for i = 1:1
    clear ph lstr;
    g = gdata(i);
    
    figure(100+i);
    clf;
    hold on;
    
    k = 1;
    if (plot_forest)
        gname = sprintf('%s/gauge%05d.txt',fc_outdir,g.id);
        if (exist(gname,'file'))
                
            tseries = importdata(gname,' ',3);

            t = tscale*(tseries.data(:,2) + tshift);
            eta = tseries.data(:,eta_idx);
            h = tseries.data(:,h_idx);
            u = vel_scale*tseries.data(:,hu_idx)./h;
            v = vel_scale*tseries.data(:,hv_idx)./h;
            speed = sqrt(u.^2 + v.^2);
            
            pvars = {eta, u, v, speed};
            pv = pvars{pvidx};
        
            ph(k) = plot(t,pv,'b.-','linewidth',2,'markersize',8);
            lstr{k} = 'ForestClaw';
            hold on;
            k = k + 1;
        else
            fprintf('File %s does not exist\n',gname);
        end
    end
                    
    if (plot_geo)
        gname_geo = sprintf('%s/gauge%05d.txt',geo_outdir,g.id);
        if (exist(gname_geo,'file'))
            tseries_geo = importdata(gname_geo,' ',3);    
            t_geo = tscale*(tseries_geo.data(:,2) + tshift);  % Shift by 10 minutes
            eta_geo = tseries_geo.data(:,eta_idx);
            h_geo = tseries_geo.data(:,h_idx);
            u_geo = vel_scale*tseries_geo.data(:,hu_idx)./h_geo;
            v_geo = vel_scale*tseries_geo.data(:,hv_idx)./h_geo;
            speed_geo = sqrt(u_geo.^2 + v_geo.^2);
            pvars_geo = {eta_geo, u_geo, v_geo, speed_geo};
            
            pv_geo = pvars_geo{pvidx};
            
            hold on;    
            ph(k) = plot(t_geo,pv_geo,'r.-','linewidth',2,'markersize',8);
            lstr{k} = 'GeoClaw';
            k = k + 1;
        else
            fprintf('File %s does not exist\n',gname_geo);
        end
    end
                        
    if (plot_obs_data)
        hold on;
        pout = plot_obs(g.id,plot_var);
        if (pout ~= 0)
            ph(k) = pout;
            lstr{k} = 'Observations';
            k = k + 1;
        end
    end
           
    xl = xlim;
    ph(k) = plot(xl,0*xl,'k');
    lstr{k} = 'Sea level';
    
    title(sprintf('Gauge %d',g.id),'fontsize',18);
    xlabel('t (hours)','fontsize',16);
    ylabel(ystr,'fontsize',16);
    set(gca,'fontsize',16);
    legend(ph,lstr);
    set(gca,'box','on');
    
    if (~plot_obs_data)
        yl = [min([pv; pv_geo]), max([pv; pv_geo])];
        ym = mean(yl);
        ys = 1.1*diff(yl)/2;
        ylim(ym + ys*[-1,1]);
    else    
        set(gca,'ylim',[-275,275]);
        % set(gca,'ylim',[-3,3]);
        set(gca,'xlim',[7.5,13]);
        % set(gca,'xtick',7:13);
        set(gcf,'position',[173   358   986   420]);
    end
       
    hold off;
    shg
       
end



end