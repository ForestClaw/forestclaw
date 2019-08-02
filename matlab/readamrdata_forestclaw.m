function [amr,t] = readamrdata_forestclaw(dim,Frame,dir)

% User-defined routine for reading in ForestClaw output.

n1 = Frame+10000;
fname = [dir, 'fort.',num2str(n1)];
fname(length(dir)+6) = 't';

if ~exist(fname,file)
    amr = {};
    t = [];
    disp(' ');
    disp(['Frame ',num2str(Frame),' (',fname,') does not exist ***']);
    disp(' ');
    return
end

% Read data from fname = 'fort.tXXXX'
fid = fopen(fname);


t = fscanf(fid,'%g',1);        fscanf(fid,'%s',1);
meqn = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
ngrids = fscanf(fid,'%d',1);   fscanf(fid,'%s',1);
fclose(fid);

% change the file name to read the q data:
fname(length(dir) + 6) = 'q';
if ~exist(fname,file)
    amr = {};
    t = [];
    return;
end

fid = fopen(fname);
disp(['Reading data from ',fname]);

for ng = 1:ngrids
    
    % read parameters for this grid:
    
    amrdata.gridno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    amrdata.level = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    amrdata.blockno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    amrdata.mpirank = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    amrdata.mx = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    amrdata.my = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    if (dim > 2)
        amrdata.mz = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    end
    
    amrdata.xlow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    amrdata.ylow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    if (dim > 2)
        amrdata.zlow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    end
    
    amrdata.dx = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    amrdata.dy = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    if (dim > 2)
        amrdata.dz = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    end
    
    % read q data:
    if (dim == 2)
        amrdata.data = fscanf(fid,'%g',[meqn,amrdata.mx*amrdata.my]);
    else
        amrdata.data = fscanf(fid,'%g',[meqn,amrdata.mx*amrdata.my*amrdata.mz]);
    end
    
    amr(ng) = amrdata;
    
end

fclose(fid);
