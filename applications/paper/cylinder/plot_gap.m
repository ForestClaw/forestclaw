function plot_gap()

showgridlines;

o = findobj('Tag','Colorbar');
delete(o);

ht = get(gca,'Title');
delete(ht);

az = -1.592298029530429e+02;
el = 2.700121742742690e+01;

zoom(3); 
view(az,el);


end