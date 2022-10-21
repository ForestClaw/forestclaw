% Axes and title
setviews;

parms = read_vars();
example = parms.example;
maxelev = parms.maxelev;
mapping = parms.mapping;



parms = read_vars();
if (mapping <= 1)
    % No mapping or Cartesian mapping
    axis([-1,1,-1,1,-1,1])
    hideslices('z',[1,3])
    daspect([1,1,1])
    view(2)
else
    % Spherical mappings with extrusion means R > 1
    s = 1 + parms.maxelev;
    axis([-1,1,-1,1,-1,1]*s)
    view(vright)
end

if qmin == qmax
    if qmax == 0
        clim([-1,1])
    else
        clim([-abs(qmax),abs(qmax)])
    end
else
    if UserVariable == 1
        % Pressure
        if mapping < 2
            clim([0.9,1.1])
        else
            clim([0.975,1.035])
        end
    else
        % Density
        clim([0.7,1.2])
        % view(2);
    end
end



fprintf("qmin = %24.16e\n",qmin);
fprintf("qmax = %24.16e\n",qmax);


daspect([1,1,1]);


if (UserVariable)
    tstr = sprintf("Presure : t = %8.2e",t);
else
    tstr = sprintf("q(%d) : t = %8.2e",mq,t);
end
% Color map and axis
colormap(parula)

title(tstr)

colorbar


% Show patch borders
showpatchborders
setpatchborderprops('linewidth',1);

shg
