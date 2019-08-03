daspect([1,1,1]);

fprintf('%15s %24.16e\n','qmin',qmin);
fprintf('%15s %24.16e\n','qmax',qmax);

if (PlotType == 1)
    colormap(parula);
    
    caxis([-1,1]);
    cv = linspace(-1,1,21);
    cv([1 end]) = [];
    % drawcontourlines(cv);
    % setcontourlineprops('linewidth',2);
    
    showpatchborders;
    if (ShowUnderOverShoots)
        u = underover;
        qlo = u.value_lower;
        qhi = u.value_upper;
        under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
        over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
        colorbar_underover(under_label,over_label);
    end
elseif (PlotType == 4)
    s = 0.1;
    axis([0,1,-1-s,1+s]);
        
    if (exist('href','var'))
        if (ishandle(href))
            delete(href);
        end
    end
    
    [amrdata_ref,tref] = readamrdata_forestclaw(2,Frame,'./fort_2nd_maxlevel6/');  
    lstyle = {'ro-','go-','bo-','mo-'};
    href = plotframe2ez(amrdata_ref,mq,lstyle,@map1d);
end


% This is used for creating vectorized PDFs
prt_tikz = false;
NoQuery = false;
% MaxFrames = 10;
if (prt_tikz)
    hidepatchborders;
    mi = 1;
    mj = 1;
    figsize = [32,32];  % Should match tikz figsize.
    maxlevel = 7;
    dpi = mi*mx*2^maxlevel/figsize(1);
    prefix = 'plot';
    caxis([-1,1]);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end



shg

