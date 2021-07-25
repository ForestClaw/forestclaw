axis([-0.5 0.5 -0.5 0.5])
view(2);
axis square

showpatchborders;
% hidepatchborders;
setpatchborderprops('linewidth',1);

fprintf('qmin = %20.16f\n',qmin);
fprintf('qmax = %20.16f\n',qmax);

colormap(parula);

%{
if (Frame == 25)
    % Use qmin, qmax from ForestClaw
    qv = [0.0702192686928361, 0.1691585563679505];
    fprintf('Caxis diff : %f\n',norm([qmin,qmax]-qv));
    caxis(qv);
else
    caxis([qmin,qmax]);
end
qv = [0.0702192686928361, 0.1691585563679505];
fprintf('Caxis diff : %f\n',norm([qmin,qmax]-qv));
caxis(qv);

% caxis([1,2]);
caxis([qmin, qmax])

dh = 1e-12;
% caxis([10, 10+dh]);
%}
qv = [0.0702192686928361, 0.1691585563679505];
caxis(qv);

set(gca,'fontsize',16);

title(sprintf('t = %g',t),'fontsize',16);

NoQuery = 0;
prt = false;
MaxFrames = 32;
if (prt)
    % hidepatchborders;
    maxlevel = 7;    %  eff. res. = 2^maxlevel
    mx = 8;
    dpi = 2^7;        % fix at 128
    figsize = (mx*2^maxlevel/dpi)*[1,1];
    prefix = 'plot_fclaw';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
end


shg

