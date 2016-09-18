function h = plot_region(ll,ur)

xp = [ll(1),ur(1),ur(1),ll(1),ll(1)];
yp = [ll(2),ll(2),ur(2),ur(2),ll(2)];
h = plot(xp,yp,'r','linewidth',3);

end
