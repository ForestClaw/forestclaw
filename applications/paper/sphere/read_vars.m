function [example,mapping,ic,kappa, period] = read_vars()

data = load('setprob.data');

example = data(1);
mapping  = data(2);
ic = data(3);
kappa = data(4);
period = data(5);


end