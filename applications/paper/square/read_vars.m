function [example,mapping,ic,alpha,center] = read_vars()

data = load('setprob.data');

example = data(1);
mapping  = data(2);
ic = data(3);
alpha = data(4);
center = data(5:6);

end