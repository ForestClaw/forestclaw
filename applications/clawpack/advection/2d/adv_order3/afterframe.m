daspect([1,1,1]);

colormap(parula);
% caxis([-1,1]);

showpatchborders;

fprintf('%15s %24.16e\n','qmin',qmin);
fprintf('%15s %24.16e\n','qmax',qmax);

if (ShowUnderOverShoots)
    u = underover;
    qlo = u.value_lower;
    qhi = u.value_upper;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    % colorbar_underover(under_label,over_label);
end
    

shg

