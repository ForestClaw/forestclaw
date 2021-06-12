if PlotType ~= 4
    cv = [0.61 : 0.02 : 1.31];
    % drawcontourlines(cv);
    
    showpatchborders
    axis([-2.5 2.5 -2.5 2.5])
    axis square
    
    setpatchborderprops('linewidth',1);
    
    set(gca,'fontsize',16);
    axis off
    
    colormap(jet);
    colorbar;    
    caxis([0.5 1.5])
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
  end;
  hold on;
  [qref,xref,p] = plotframe1ez(rad1d,mq,'b-');
  set(p,'LineWidth',2,'markersize',10);

  [h2d,leg_str] = getlegendinfo(0);
  legend([h2d,p],{leg_str{:},'1d reference soln'});
  hold off;

end


shg
