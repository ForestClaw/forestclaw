% Axes and title
setviews;

fprintf("qmin = %24.16e\n",qmin);
fprintf("qmax = %24.16e\n",qmax);

% Application parameters (written in overpressure_user.cpp)
parms = read_vars();
example = parms.example;
maxelev = parms.maxelev;
mapping = parms.mapping;


% Set axis
parms = read_vars();
switch mapping
    case 0
        % No mapping or Cartesian mapping    
        axis([-1,1,-1,1,-1,1])
        daspect([1,1,1])
        view(3)
    case 1
        % Cartesian mapping    
        axis([-1,1,-1,1,parms.minz,parms.maxz])
        daspect([1,1,1])
        view(3)
    case 2
        % Spherical mappings with extrusion means R > 1
        s = 1 + parms.maxelev;
        axis([0,1,-1,1,-1,1]*s)
        view(3)
    case {3,4}
        % Spherical mappings with extrusion means R > 1
        s = 1 + parms.maxelev;
        axis([-1,1,-1,1,-1,1]*s)
        view(3)
    otherwise
        error("No mapping specified.")
end
set(gca,'clipping','off')
set(gca,'box','on')



% Set color axis
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
            clim([0.8,1.2])
        else
            clim([0.975,1.125])
        end
    else
        % Density
        clim([-1,1]*1e-3)
    end
end

if (UserVariable)
    tstr = sprintf("Pressure : t = %8.2e",t);
else
    tstr = sprintf("q(%d) : t = %8.2e",mq,t);
end
title(tstr);

% Color map and axis
colormap(parula)
colorbar


% cv = linspace(0.8,1.2,24);
% drawcontourlines(cv);
showpatchborders
    

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
