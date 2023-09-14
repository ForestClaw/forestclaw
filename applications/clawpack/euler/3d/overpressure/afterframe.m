% Axes and title
setviews;

fprintf("qmin = %24.16e\n",qmin);
fprintf("qmax = %24.16e\n",qmax);

% Application parameters (written in overpressure_user.cpp)
parms = read_vars();
example = parms.example;
maxelev = parms.maxelev;
mapping = parms.mapping;

% -------------------------------------
% Set axis for mapping
% -------------------------------------
parms = read_vars();
switch mapping
    case 0
        % No mapping
        axis([-1,1,-1,1,-1,1])
        daspect([1,1,1])
        view(3)
    case 1
        % Cartesian mapping    
        if parms.scale_bump > 0
            axis([-2,2,-2,2,parms.minz,parms.maxz])
        else
            axis([-1,1,-1,1,parms.minz,parms.maxz])
        end
        daspect([1,1,1])
        view(3)
    case 2
        % Spherical mappings with extrusion means R > 1
        s = 1 + parms.maxelev;
%         axis([-1,1,-1,1,-1,1]*s)
        axis([0,2.5,-2,2,-2,2])    
    case {3,4}
        % Spherical mappings with extrusion means R > 1
        s = 1 + parms.maxelev;
        axis([-1,1,-1,1,-1,1]*s)
        set(gcf,'clipping','off')
        view(vright)
%         axis([-1.25,1.25,0,1.25,-1.25,1.25])        
%         set(gca,'clipping','on')
%         view([42.8,12.03])
    otherwise
        error("No mapping specified.")
end
set(gca,'box','on')

% showgridlines
showpatchborders
    
% hideslices
showslices('z',3)

set(gca,'clipping','off')

% -------------------------------------
% Color axis and color map
% -------------------------------------
if qmin == qmax
    if qmax == 0
        clim([-1,1])
    else
         clim([-abs(qmax),abs(qmax)])
    end
else
    if UserVariable == 1
        if strcmp(UserVariableFile,'pressure') == 1
            % Pressure
            clim([0.9,1.1])
        end
    else
        clim([0.95,1.05])
    end
end
% Color map and axis
colormap(parula)
colorbar

cl = clim;
cv = linspace(cl(1),cl(2),22);    % Try to skip 1
drawcontourlines(cv)
setcontourlineprops('linewidth',1);


% -------------------------------------
% Title and labels
% -------------------------------------

if (UserVariable)
    tstr = sprintf("Pressure : t = %8.2e",t);
else
    tstr = sprintf("q(%d) : t = %8.2e",mq,t);
end
title(tstr);

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
