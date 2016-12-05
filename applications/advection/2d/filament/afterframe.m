s = 0.05;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);
caxis([0 1]);

N = 500;
if (t > 0)
    [xout,yout] = filament_soln(N,t);
else
    th = linspace(0,2*pi,N+1);
    xout = 0.25*cos(th) + 0.5;
    yout = 0.25*sin(th) + 1;
end
hold on;
plot(xout,yout,'k','linewidth',2);
hold off;

showpatchborders(1:7);
setpatchborderprops('linewidth',1)

view(2);
axis off
shg;


NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 41;
  filename = sprintf('filament_%04d.png',Frame)
  print('-dpng',filename);
end

NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 32;
  filename = sprintf('filament%03d.png',Frame);
  fprintf('Printing %s ...\n',filename);
  print('-dpng',filename);
end

%%
prt = false;
NoQuery = 0;
if (prt)
    MaxFrames = 41;
    axis([0 2 0 2]);
    axis off;
    delete(get(gca,'title'));
    hidepatchborders;
    figsize = [4,4];  % Should match size set in options
%     set(gcf,'papersize',figsize);
%     set(gca,'position',[0 0 1 1]);
%     set(gcf,'paperposition',[0 0 figsize]);

    % Use this with 'export_fig'
     set(gca,'position',[0 0 1 1]);
     set(gcf,'units','inches');
     set(gcf,'position',[1 7 figsize]);

    % Start printing
    id = input('Input id to use : ');
    if (~isempty(id) | id == 999)

        % fname_soln_tikz = sprintf('filament_tikz_%04d.tex',Frame);
        % fname_soln_tikz_dir = sprintf('results_%03d/filament_tikz_%04d.tex',id,Frame);
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
        fname_png = sprintf('results_%03d/%s_%04d.png',id,fname_prefix,Frame);
        if (exist(fname_png))
            str = sprintf('Overwrite file %s (y/[n]) ? ',fname_png);
            yn = input(str,'s');
            if (isempty(yn))
                yn = 'n';
            end
        end

        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname_png);
            % print('-dpng','-r256',fname_png);
            export_fig('-dpng','-transparent','-r256',...
               '-a1','-p0','-nocrop',fname_png);
            amrclaw = 0;
            % Frame,fname_prefix,id,amrclaw,extra_file
            create_tikz_plot(Frame,fname_prefix,id,amrclaw);
        end

    end

end

clear afterframe
clear mapc2m
