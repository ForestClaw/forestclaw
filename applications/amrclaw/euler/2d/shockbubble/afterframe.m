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
    set(gcf,'papersize',[4,1]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 4 1]);
    id = input('Input id to use : ');
    if ~isempty(id)
    
        % With mesh
        setpatchborderprops('linewidth',1);
        hidegridlines;
        if (PlotType == 3)
            fname = sprintf('results_%03d/fc_sb_schlrn_%02d.png',id,Frame);            
        else
            fname = sprintf('results_%03d/fc_sb_%02d.png',id,Frame);
        end
        yn = 'y';
        if (exist(fname))
            str = sprintf('Are you sure you want to overwrite %s (y/[n]) ? ',...
                fname);            
            yn = input(str,'s');               
            if isempty(yn)
                yn = 'n';
            end
        end      
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname);
            print('-r1024','-dpng',fname);
        end
        
        % No mesh
        hidegridlines;
        hidepatchborders;
        if (PlotType == 3)
            fname = sprintf('results_%03d/fc_sb_schlrn_nomesh_%02d.png',id,Frame);            
        else
            fname = sprintf('results_%03d/fc_sb_nomesh_%02d.png',id,Frame);
        end
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname);
            print('-r1024','-dpng',fname);
        end
                
    end
    
end

shg;

clear afterframe
clear mapc2m
clear parallelpartitions