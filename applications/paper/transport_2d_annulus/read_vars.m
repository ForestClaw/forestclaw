function [example,A,rinit,init_location, beta,theta,vcart, ...
    freq,cart_speed] = read_vars()


d = load('_output/mapping.dat');
example = d(1);
A = d(2);
rinit = d(3);
init_location = d(4:5);
beta = d(6);
theta = d(7:8);
vcart = d(9:10);
cart_speed = d(11);
freq = d(12);

end