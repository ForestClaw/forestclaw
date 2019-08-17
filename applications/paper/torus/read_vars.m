function [example,ic,alpha, beta,revs_per_sec] = read_vars()

data = load('setprob.data');

example = data(1);
ic = data(2);
alpha = data(3);
beta = data(4);
revs_per_sec = data(5);


end