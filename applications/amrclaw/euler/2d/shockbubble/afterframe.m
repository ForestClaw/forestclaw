if (PlotType ~= 3)
    colormap(jet)
end


axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4e\n','qmin',qmin);
fprintf('%10s : %12.4e\n','qmax',qmax);

if (PlotType == 3)
    caxis([0 200]);
else
    caxis([0.1 2.81]);
end

showpatchborders
setpatchborderprops('linewidth',1);
% showgridlines(1:5);
% hidepatchborders;
set(gca,'fontsize',16);

tstr = sprintf('ForestClaw : t = %12.4f',t);
title(tstr,'fontsize',16);

if (PlotParallelPartitions==1)
    showpatchborders;
end

prt = true;
NoQuery = 0;
if (prt)
    MaxFrames = 31;
    axis off;
    delete(get(gca,'title'));
    figsize = [8,2];  % Should match size set in options
    set(gcf,'papersize',figsize);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 figsize]);  
    
    % Start printing
    id = input('Input id to use : ');
    if ~isempty(id)
    
        % No mesh
        hidegridlines;
        hidepatchborders;
        if (PlotType == 3)
            fname_prefix = sprintf('fc_sb_schlrn',Frame);            
        else
            fname_prefix = sprintf('fc_sb',Frame);
        end
        yn = 'y';
        fname_png = sprintf('results_%03d/%s_%04d.png',id,fname_prefix,Frame);
        if (exist(fname_png))
            str = sprintf('Overwrite file %s (y/n) ? ',fname_png);
            yn = input(str,'s');
            if (isempty(yn))
                yn = 'n';
            end
        end
                
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname_png);
            print('-r512','-dpng',fname_png);
            plotgrid = 1;
            create_tikz_plot(id,Frame,fname_prefix,1);
        end
                
    end
    
end

shg;

clear afterframe
clear mapc2m
clear parallelpartitions