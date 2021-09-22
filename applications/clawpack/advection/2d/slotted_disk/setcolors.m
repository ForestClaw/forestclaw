function setcolors(p,x,y,z,q)

global ash_cm_limit ash_cm_N;

ash_cm_limit = [4, 100];
ash_cm_N = 4;

N = ash_cm_N;
qv = ash_cm_limit;

b = [102, 178,255]/255;
c1 = b;
d = linspace(0,1,N+1)';
c2 = kron(d(1:end-1),[1,1,1]);
c3 = [1,1,1];
cm = [c1; c2; c3];

idx = zeros(size(q));

m1 = q < qv(1);
idx(m1) = 1;
for i = 1:N
    q1 = qv(1) + d(i)*(qv(2)-qv(1));
    q2 = qv(1) + d(i+1)*(qv(2)-qv(1));
    m = q1 <= q & q <= q2;
    idx(m) = i+1;
end
m3 = q > qv(2);
idx(m3) = N+2;

colormap(cm);

set(p,'cdata',idx);
fv_idx = get(p,'FaceVertexCData');
set(p,'FaceVertexCData',fv_idx);

set(p,'cdatamapping','direct');
set(p,'FaceColor','flat');


end