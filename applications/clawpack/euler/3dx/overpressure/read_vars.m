function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.mapping = data(2);
parms.manifold = data(3);
% ...
parms.longitude = data(15:16);
parms.latitude = data(17:18);
parms.maxelev = data(19);
parms.minz = data(20);
parms.maxz = data(21);
parms.midz = data(22);
parms.scale_bump = data(23);

parms.scale = data(24:25);
parms.shift = data(26:27);



end