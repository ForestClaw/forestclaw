function [example,A,rinit,beta,theta,freq,cart_speed] = read_vars();


d = load('mapping.dat');
example = d(1);
A = d(2);
rinit = d(3);
beta = d(4);
theta = d(5:6);
freq = d(7);
cart_speed = d(8);

end