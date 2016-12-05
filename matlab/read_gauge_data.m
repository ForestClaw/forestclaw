function gdata = read_gauge_data(gnum)

fname = sprintf('gauge%05d.txt',gnum);
if (~exist(fname))
    error(sprintf('Gauge file %s does not exist.',fname));
end
fid = fopen(fname,'r');
k = 1;
while (true)
    l = fgetl(fid);
    if (~ischar(l))
        break
    elseif (~strcmp(l(1),'#'))
        gdata(k,:) = sscanf(l,'%d %e %e %e %e %e',Inf);
        k = k + 1;
    end
end


end