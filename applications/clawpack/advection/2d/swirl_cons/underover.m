function uo = underover()
%UNDEROVER returns structure for assigning colors to under/overshoots
%
%    UO = UNDEROVER.  To have special highlighting for under/over shoot
%    values, the user should first set the ShowUnderOverShoots flag in
%    SETPLOT2 or the SETPLOT3 files.
%    Then, make a local copy of this file and supply values for the the
%    following fields in the UO structure :
%
%             uo.color_under     the rgb triple used to color undershoots
%             uo.color_over      the rgb triple used to color overshoots
%             uo.color_nan       the color to use for NAN values
%             uo.value_lower     lower value of exact solution
%             uo.value_upper     upper value of exact solution
%             uo.tol             tolerance for under/over shoots
%             uo.colormap        colormap for values in [qlow-tol,qhi+tol];
%
%    Values are considered in the overshoot/undershoot region if they satisfy
%
%                    q < value_lower - tol      (undershoot)
%    or
%                    q > value_upper + tol      (overshoot)
%
%    Values not in the under/overshoot region are linearly scaled into the
%    colormap.  Values that lie in [qlow-tol,qlow] are clamped to qlow, and
%    values in [qhi, qhi + tol] are clamped to qhi.
%
%    The under/over shoot colorbar only takes affect if you have set
%    'ShowUnderOverShoots' in SETPLOT2.m.
%
%    Example :  First, create underover structure :
%    function uo = underover()
%
%    colormap('default');
%    cm = colormap;
%
%    uo = struct('color_under',[0 1 1],...  % cyan
%                'color_over',[1 0 1],...   % magenta
%                'color_nan',[1 1 1],...    % white
%                'value_lower',0, ...       % theoretical minimum
%    	         'value_upper',1,...        % theoretical maximum
%    	         'tol',1e-4,...             % acceptable numerical over/undershoot
%    	         'colormap',cm);            % Colors for non under/over shoot values
%
%
%    Example : To make sure that your under/over shoot color map is used :
%    % In local copy of 'setplot2.m'
%    ShowUnderOverShoots = 1;
%
%    % In local copy of 'afterframe.m'
%    qlo = 0;
%    qhi = 1;
%    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
%    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
%    colorbar_underover(under_label,over_label);
%
%    See also SETCOLORS, COLORBAR_UNDEROVER.
%

cm = yrbcolormap;

uo = struct('color_under',[0 1 1],...  % cyan
    'color_over',[1 0 1],...   % magenta
    'color_nan',[1 1 1],...    % white
    'value_lower',0, ...       % theoretical minimum
    'value_upper',2.738465249609046,...   % theoretical maximum
    'tol',eps(1),...             % acceptable numerical over/undershoot
    'colormap',cm);            % Colors for non under/over shoot values


