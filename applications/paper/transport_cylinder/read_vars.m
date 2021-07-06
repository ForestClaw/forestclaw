function [example,ic,refine_pattern, R, H,xc0,yc0,r0, ...
    revs_per_sec,v_speed] = read_vars()

if (exist('_output/mapping.dat','file'))
    data = load('_output/mapping.dat');
else
    error('File ''_output/mapping.dat'' not found.');
end

example = data(1);
ic = data(2);
refine_pattern = data(3);
R = data(4);
H = data(5);
xc0 = data(6);
yc0 = data(7);
r0 = data(8);
revs_per_sec = data(7);
v_speed = data(8);


end