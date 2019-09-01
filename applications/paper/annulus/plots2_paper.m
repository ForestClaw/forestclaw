function plots2_paper()

[example,A,rinit,beta,theta,freq,cart_speed,refine_pattern] = read_vars();

axis([-0.1, 0.1, 0.6, 0.76])
showgridlines;
showpatchborders;
setpatchborderprops('linewidth',2);

set(gca,'fontsize',16);
th = get(gca,'title');
set(th,'fontsize',18);

if refine_pattern == 0
    print('-dpng','constant_theta.png');
elseif refine_pattern == 1
    print('-dpng','constant_r.png');
end

end