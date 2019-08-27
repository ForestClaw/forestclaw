function plot_gap()

[example,ic,refine_pattern, alpha, beta,init_radius, ...
    revs_per_sec,cart_speed,theta_range,phi_range] = read_vars();

showgridlines
showpatchborders
setpatchborderprops('linewidth',2)

setslicecolor([1,1,1]*0.8)
o = findobj('type','colorbar');
if ishandle(o)
    delete(o);
end

axis off
delete(get(gca,'title'))

set(gca,'Clipping','off');
% set(gca,'dataaspectratiomode','auto')
daspect([1,1,1]);

v = [35.70, 19.54];
view(v)



end