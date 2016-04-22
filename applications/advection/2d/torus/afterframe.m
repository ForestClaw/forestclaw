setviews;

global map isflat;

alpha = 0.4;
s = 1e-2;    
if strcmp(map,'nomap') || strcmp(map,'brick')
    % axis([-1 1 -1 1]);
    axis([0 1 0 1])
    % axis image;
else
    if strcmp(map,'cart')
        alim = [-1 1];
    elseif strcmp(map,'latlong')
        alim = [-1 1];
    else
        alim = [-1-alpha,1+alpha];
    end
    alim = alim + [-s s];
    axis([alim alim]);
    daspect([1 1 1]);
    view(vtop)
end

if PlotParallelPartitions ~= 1
    yrbcolormap;
end
showpatchborders(1:10);
caxis([0,1])
qlo = 0;
qhi = 1;
under_label = sprintf('0 - %7.1e',qlo-qmin);
over_label = sprintf('1 + %7.1e',qmax-qhi);
fprintf('%6s %12s\n','qmin',under_label);
fprintf('%6s %12s\n\n','qmax',over_label);

if (ShowUnderOverShoots)
    qlo = 0;
    qhi = 1;
    colorbar_underover(under_label,over_label);
end
daspect([1,1,1]);

if (isflat)
    view(2);
else
    hidepatchborders;
    view(3);
    % camlight;
    showpatchborders;
    hold on;
    R = 0.5;  % revs per second
    alpha = 0.4;
    period = 16;
    [xp,yp,zp] = torus_soln(t,alpha,R,period);    % camlight;
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
end

setpatchborderprops('linewidth',1)
hidepatchborders(9)

%%
cm = ...
    [166	206	227;
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
    177	89	40]/255;    


%%
NoQuery = 0;
prt = true;
if (prt)
  MaxFrames = 128;
  axis([0 1 0 1]);
  axis off;
  figsize = [4,4];  % Should match size set in options
  delete(get(gca,'title'));
  fname_png = sprintf('torus_pp_%04d.png',Frame);
  % print('-dpng',filename);
  set(gca,'position',[0 0 1 1]);
  set(gcf,'units','inches');
  figsize = [4,4];  % Should match size set in options
  set(gcf,'position',[1 7 figsize]);
  export_fig('-dpng','-r256',...
      '-a1','-p0','-nocrop',fname_png);
end

%%
prt = false;
% NoQuery = 0;
if (prt)
    MaxFrames = 41;
    axis([0 1 0 1]);
    axis off;
    delete(get(gca,'title'));
%     hidepatchborders;
    figsize = [4,4];  % Should match size set in options
%     set(gcf,'papersize',figsize);
%     set(gca,'position',[0 0 1 1]);
%     set(gcf,'paperposition',[0 0 figsize]);

    % Use this with 'export_fig'
     set(gca,'position',[0 0 1 1]);
     set(gcf,'units','inches');
     set(gcf,'position',[1 7 figsize]);

    % Start printing
%     fname_soln_tikz = sprintf('torus_tikz_%04d.tex',Frame);
%     fname_soln_tikz_dir = sprintf('torus_tikz_%04d.tex',Frame);
    % create_filament_soln_tikz(fname_soln_tikz_dir,xout,yout,figsize,1024,1024);


    % No mesh
    hidegridlines;
    % hidepatchborders;
    if (PlotType == 3)
        fname_prefix = sprintf('fc_adv_schlrn',Frame);
    else
        fname_prefix = sprintf('fc_adv',Frame);
    end
    yn = 'y';
    fname_png = sprintf('%s_%04d.png',fname_prefix,Frame);
    if (exist(fname_png))
        str = sprintf('Overwrite file %s (y/[n]) ? ',fname_png);
        yn = input(str,'s');
        if (isempty(yn))
            yn = 'n';
        end
    end

    if (strcmp(lower(yn),'y') == 1)
        fprintf('Printing %s\n',fname_png);
        % print('-dpng','-r1024',fname_png);
        export_fig('-dpng','-transparent','-r1024',...
            '-a1','-p0','-nocrop','-native',fname_png);
        amrclaw = 0;
        create_tikz_plot(Frame,fname_prefix,amrclaw,fname_soln_tikz);
    end

end


shg

clear afterframe;
clear mapc2m;
clear mapc2m_torus;
clear torus_soln;
clear parallelpartitions;
