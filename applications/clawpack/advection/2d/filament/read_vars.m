function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.alpha = data(2);
parms.center = data(3:4);

end