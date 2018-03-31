function plot_tikz_fig(Frame,figsize,prefix,dpi)
% PLOT_TIKZ_FIG prints figure for use with tikz
%
% PLOT_TIKZ_FIG(Frame,figsize,prefix) prints a figure with name
% '<prefix>.NNNN.png', where NNNN is the zero-padded frame number,
% and 'prefix' is a string.   The figsize is used to size the figure
% appropriately for a tikz mesh created with Latex.
%
% PLOT_TIKZ_FIG(Frame,figsize,prefix,dpi) allows for an optional 
% resolution DPI (dots-per-inch).  To match tikz figure, set res=mx*2^maxlevel.
%
% This is to be used with the ForestClaw tikz plotting option.
%
% See also clawgraphics.

    dpi_test = dpi;    
    scaling_factor = 1;
    while dpi_test > 128
        dpi_test = dpi_test/2;
        scaling_factor = 2*scaling_factor;
    end        
    if (scaling_factor > 1)
        fprintf('\n');
        fprintf('----> Warning : Specified dpi > 128; dpi = %d.  Consider reducing\n', dpi);
        fprintf('      the dpi and increasing figure size by a factor of %d.\n',scaling_factor);
        fprintf('\n');
    end

    fs = figsize(:)';
    axis off
    hidepatchborders;
    delete(get(gca,'title'));
    set(gcf,'papersize',fs);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'paperposition',[0 0 fs]);
    fname = sprintf('%s_%04d.png',prefix,Frame);
    
    % Match print resolution to computational resolution
    if (nargin == 4)        
        resstr = sprintf('-r%d',dpi);
        print('-dpng',resstr,fname);  
    else
        print('-dpng',fname);
    end
end
