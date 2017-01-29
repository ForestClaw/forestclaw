if (PlotType ~= 3)
    colormap(jet)
end

%%
axis([0 1.5 -0.75 0.75])
axis equal;
axis tight
axis off

%%
fprintf('%10s : %12.4e\n','qmin',qmin);
fprintf('%10s : %12.4e\n','qmax',qmax);

if (PlotType == 3)
    caxis([0 200]);
elseif mq == 5
    % 0=y; 1=b; g=2; r=3;
    b = [0 0 1];
    g = [0 0.8 0];
    r = [1 0 0];
    y = [1 1 0];
    s = linspace(0,1,11)';
    yb = kron(y,1-s) + kron(b,s);
    bg = kron(b,1-s) + kron(g,s);
    gr = kron(g,1-s) + kron(r,s);
    ry = kron(r,1-s) + kron(y,s);
    cm = [y; yb; b; bg; g; gr; r];
    cm = [r; ry; y; yb; b; bg; g];
    colormap((cm));
    
    
    caxis([0 3])
    % yrbcolormap;
else
    colormap(jet)
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

showpatchborders;
% colormap(white)

prt = false;
NoQuery = 0;
if (prt)
    MaxFrames = 1000;
    axis off;
    delete(get(gca,'title'));
    figsize = [8,8];  % Should match size set in options
    set(gcf,'papersize',figsize);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 figsize]);

    % Use this with 'export_fig'
%     set(gca,'position',[0 0 1 1]);
%     set(gcf,'units','inches');
%     set(gcf,'position',[0 0 figsize]);

    % Start printing
    % No mesh
    hidegridlines;
    hidepatchborders;
    if (PlotType == 3)
        fname_prefix = sprintf('fc_sb_schlrn_%04d.png',Frame);
    else
        fname_prefix = sprintf('fc_sb',Frame);
    end
    yn = 'y';
    id = input('Input id : ');   
    if (isempty(id))
        id = 999;
    end
    fname_png = sprintf('results_%03d/%s_%04d.png',id,fname_prefix,Frame);
    if (exist(fname_png,'file'))
        str = sprintf('Overwrite file %s (y/[n]) ? ',fname_png);
        yn = input(str,'s');
        if (isempty(yn))
            yn = 'n';
        end
    end
    if (strcmp(lower(yn),'y') == 1)
        fprintf('Printing %s\n',fname_png);
        print('-dpng','-r256',fname_png);
        %             export_fig('-dpng','-transparent','-r512',...
        %                 '-a1','-nocrop',fname_png);
        create_tikz_plot(Frame,fname_prefix,id);
    end
end

shg;

clear afterframe
clear mapc2m
clear parallelpartitions
