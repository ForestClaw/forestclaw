function create_filament_soln_tikz(fname,xout,yout,figsize,mxf,myf)

ax = 0;
bx = 2;
ay = 0;
by = 2;

lx = (bx-ax);
ly = (by-ay);

dxf = lx/mxf;
dyf = ly/myf;

sx = figsize(1)/mxf;
sy = figsize(2)/myf;

fid = fopen(fname,'w');
fprintf(fid,'\\begin{tikzpicture}[x=%18.16fin,y=%18.16fin]\n',sx,sy);
fprintf(fid,'\\useasboundingbox (0,0) rectangle (%d,%d);\n',mxf,myf);



for i = 1:length(xout)-1,
    x1 = (xout(i)-ax)/dxf;
    y1 = (yout(i)-ay)/dyf;
    x2 = (xout(i+1)-ax)/dxf;
    y2 = (yout(i+1)-ay)/dyf;
    fprintf(fid,'\\draw[thick,black] (%8.2f, %8.2f) -- (%8.2f, %8.2f);\n',x1,y1,x2,y2);
end
fprintf(fid,'\\end{tikzpicture}\n');