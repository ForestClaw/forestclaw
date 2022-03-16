setviews;

% yrbcolormap;
colormap(parula);

fprintf('qmin = %16.8e\n',qmin);
fprintf('qmax = %16.8e\n',qmax);


if (UserColorMapping > 0)
    % showgridlines(1:4);
    showpatchborders;
    setpatchborderprops('linewidth',1);
    hidepatchborders(7);
    view([65.5,12]);
%{
    global ash_cm_limit ash_cm_N;
    qv = ash_cm_limit;
    N = ash_cm_N;

    o = findobj('Tag','Colorbar');
    clabel = cell(N+3,1);
    clabel{1} = 0;
    clabel{2} = qv(1);
    clabel{end-1} = qv(2);
    if (qmax <= qv(2))
        clabel{end} =  '';
    else
        clabel{end} = sprintf('%.0f',qmax);
    end
    set(o,'yticklabel',clabel);
    set(o,'ticklength',[0])
    set(o,'fontsize',16,'fontweight','bold')
%}
else
    showpatchborders;
    setpatchborderprops('linewidth',1);
    caxis([-0.5,1]);
end

colorbar

axis image;
set(gca,'clipping','off')
daspect([1 1 1]);
axis off;

view([117.50,4]);

delete(get(gca,'title'));


NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 10;  
  filename = framename(Frame,'sdisk000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng','-r600',filename);
end

shg
