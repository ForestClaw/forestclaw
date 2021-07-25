if PlotType ~= 4
    cv = [0.61 : 0.02 : 1.31];
    % drawcontourlines(cv);
    
    showpatchborders
    axis([-2.5 2.5 -2.5 2.5])
    axis square
    
    setpatchborderprops('linewidth',1);
    
    set(gca,'fontsize',16);
    axis off

    if (mq > 1)
        caxis([-1,1])
    else
        caxis([0.5 1.5])
    end
    
    colormap(jet);
    colorbar;    
    view(2)
end

if PlotType==4
  axis([0 2.5 0 2.1])
  [rad1d,tref] = readamrdata(1,Frame,'./1drad/');
  if isempty(rad1d)
    disp('Run xclaw in rad1d directory to make reference solution')
    return
  end
  if (abs(tref-t) > 1e-8)
    error('times are not compatible');
  end
  hold on;
  [qref,xref,p] = plotframe1ez(rad1d,mq,'k-');
  set(p,'LineWidth',1.5,'markersize',10);

  [h2d,leg_str] = getlegendinfo(0);
%   lh = legend([h2d,p],{leg_str{:},'1d reference soln'});
%   set(lh,'AutoUpdate','off');  
  hold off
  
  hold on
  xm = (xref(1:end-1) + xref(2:end))/2;
  pgrad = plot(xm(:),diff(qref(:))./diff(xref(:)),'r','linewidth',2);
  axis([0,2.5,-2,6]);
  
  lh = legend([h2d,p,pgrad],{leg_str{:},'1d reference soln','Gradient'});
  set(lh,'AutoUpdate','off');  
    
  hold off;

end


shg
