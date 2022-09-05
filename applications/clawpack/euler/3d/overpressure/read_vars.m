function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.mapping = data(2);
parms.manifold = data(3);
% ...
parms.longitude = data(15:16);
parms.latitude = data(17:18);
parms.maxelev = data(19);

end