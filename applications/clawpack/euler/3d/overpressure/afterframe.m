% Axes and title

parms = read_vars();
example = parms.example;
maxelev = parms.maxelev;
mapping = parms.mapping;




parms = read_vars();
if (mapping <= 1)
    axis([-1,1,-1,1,-1,1])
    hideslices('z',[1,3])
    daspect([1,1,1])
    view(3)
else
    % Spherical mappings with extrusion means R > 1
    s = 1 + parms.maxelev;
    axis([-1,1,-1,1,-1,1]*s)
    view([119,15])
end

if qmin == qmax
    if qmax == 0
        caxis([-1,1])
    else
        caxis([-abs(qmax),abs(qmax)])
    end
else
    if UserVariable == 1
        % Pressure
        caxis([0.95,1.15])
    else
        % Density
        caxis([0.5,1.3])
        % caxis([-1,1]*1e-3);  % For errors
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
colorbar

title(tstr)


% Show patch borders
showpatchborders
setpatchborderprops('linewidth',1);
