% Axes and title
setviews;

parms = read_vars();
example = parms.example;
maxelev = parms.maxelev;
mapping = parms.mapping;



parms = read_vars();
if (mapping <= 1)
    % No mapping or Cartesian mapping
    
    axis([-1,1,-1,1,parms.minz,parms.maxz])
    hideslices('z')
    daspect([1,1,1])
    view(vfront)
    % view(3)
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
        clim([-1,1]*1e-3)
        showslices('z');
        view(2);
    end
end

fprintf("qmin = %24.16e\n",qmin);
fprintf("qmax = %24.16e\n",qmax);


daspect([1,1,1]);


if (UserVariable)
    tstr = sprintf("Pressure : t = %8.2e",t);
else
    tstr = sprintf("q(%d) : t = %8.2e",mq,t);
end
title(tstr);

% Color map and axis
colormap(parula)
colorbar


compare_files = true;
if compare_files
    clim([-1,1]*1e-4)
    showslices('z');
    view(2);
else
    cv = linspace(0.8,1.2,24);
    drawcontourlines(cv);
    % Show patch borders
    showpatchborders
end
    




if mapping <= 1
    prt = false;
    if prt
        setpatchborderprops('linewidth',2);
        hideslices('x')
        fstr = sprintf("bump_%02d.png",Frame);
        fprintf("Printing file '%s'\n",fstr);
        print(fstr,'-dpng');
    end
end


shg
