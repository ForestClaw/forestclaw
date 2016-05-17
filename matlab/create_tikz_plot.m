function create_tikz_plot(id,Frame,fname_prefix,amrclaw,extra_file)

% \documentclass{standalone}
%
% \usepackage{tikz}
%
% \newcommand{\plotgrid}[1]{}
%
% % Needed only for AMRClaw
% \newcounter{level}
% \makeatletter
% \renewcommand\thelevel{\two@digits{\value{level}}}
% \makeatother
% \newcommand{\addlevel}{\stepcounter{level}\includegraphics{amr_sb_0000_\thelevel.png}}
%
% \begin{document}
% \input{tikz.0000.tex}
% \end{document}

if nargin < 5
    extra_file = [];
    if (nargin < 4)
        amrclaw = 0;        
    end
end

for p = 0:1,
    if p == 0
        % Basic plot with no mesh
        fname_tex = sprintf('%s_%04d.tex',fname_prefix,Frame);
    else
        fname_tex = sprintf('%s_mesh_%04d.tex',fname_prefix,Frame);
    end
    fprintf('Printing %s\n',fname_tex);
    fid = fopen(fname_tex,'w');
    fprintf(fid,'\\documentclass{standalone}\n');
    fprintf(fid,'\\usepackage{tikz}\n');
    if p == 0
        fprintf(fid,'\\newcommand{\\plotgrid}[1]{}\n');
    else
        fprintf(fid,'\\newcommand{\\plotgrid}[1]{#1}\n');
    end
    if (amrclaw == 1)
        % We need to include a PNG file for each level, and layer the levels,
        % and grid plotting in Tikz. The extra latex commands create file
        % names with a zero-padded level number, e.g. amr_sb_0000_01.png.
        % Each call to \addlevel increments the counter 'level'.
        fname_png = sprintf('%s_%04d_\\thelevel.png',fname_prefix,Frame);
        fprintf(fid,'\\newcounter{level}\n');
        fprintf(fid,'\\makeatletter\n');
        fprintf(fid,'\\renewcommand\\thelevel{\\two@digits{\\value{level}}}\n');
        fprintf(fid,'\\makeatother\n');
        fprintf(fid,'\\newcommand{\\addlevel}{\\stepcounter{level}\\includegraphics{%s}}\n',fname_png);
    else
        fname_png = sprintf('%s_%04d.png',fname_prefix,Frame);
        fprintf(fid,'\\newcommand{\\figname}{%s}\n',fname_png);
    end
    
    
% \begin{tikzpicture}
% \node (mesh) {\input{tikz.0000.tex}};
% \node (mesh) {\input{filament_tikz_0000.tex}};
% \end{tikzpicture}
    
    fprintf(fid,'\\begin{document}\n');
    fprintf(fid,'\\begin{tikzpicture}\n');
    fprintf(fid,'\\node (mesh) {\\input{tikz.%04d.tex}};\n',Frame);
    %fprintf(fid,'\\input{tikz.%04d.tex}\n',Frame);
    efile = sprintf('%s',extra_file);
    if (exist(efile,'file'))
        fprintf(fid,'\\node (mesh) {\\input{%s}};\n',extra_file);
        % fprintf(fid,'\\input{%s}\n',extra_file);
    end
    fprintf(fid,'\\end{tikzpicture};\n');
    fprintf(fid,'\\end{document}\n');
    fclose(fid);
end