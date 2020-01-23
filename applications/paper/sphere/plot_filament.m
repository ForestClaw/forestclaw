function plot_filament(dir)

close all

global show_degs

% Limiters
% 0 - no limiter
% 1 - Minmod
% 2 - Superbee
% 3 - Van Leer
% 4 - Montonized Centered


nolim = false;
withlim = [false, false, false, true];  % Use limiter ID
levels = [true, true, true, true];      % Levels 1,2,3,4
uniform = false;
adapt = true;

info = struct('nolim',nolim,'withlim',withlim,'levels',levels,...
    'uniform',uniform,'adapt',adapt);

show_degs = true;
plot_last = false;

lsu = ':';  % line style for uniform
lsa = '-'; % line style for adaptive

mx = 32;
nf = 1;    % file counter (incremented each time a file is added)

nfiles_max = 30;   % Maximum number of files
% ------------------------- Initial conditions ----------------------------
fname = 'level2_uniform_nolim_t0.dat';
lstr = sprintf('t=0.0');
s(nf) = add_to_plot(dir,fname,lstr,nf,'r','.',lsu,30,false);
nf = nf+1;

% -------------------------------- Test -----------------------------------
add2plot = false;
fname = 'filament.dat';
lstr = sprintf('TEST');
s(nf) = add_to_plot(dir,fname,lstr,nf,'k','.','-',30,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL = 1 --------------------------------
maxlevel = 1;
sym_color = 'r';

lim = 0;
add2plot = nolim & levels(maxlevel) & info.uniform;
fname = 'level1_uniform_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

lim = 1;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level1_uniform_withlim_1_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

lim = 2;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level1_uniform_withlim_2_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

lim = 3;
add2plot = withlim(lim) & levels(maxlevel)  & info.uniform; 
fname = 'level1_uniform_withlim_3_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level1_uniform_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ----------------

lim = 0;
add2plot = nolim & levels(maxlevel) & info.adapt;
fname = 'level1_adapt_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,15,add2plot);
nf = nf+1;

% Adaptive runs only done for limiter=4
lim =4;
add2plot = withlim(4) & levels(maxlevel) & info.adapt;
fname = 'level1_adapt_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,8,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=2 ----------------------------------
maxlevel=2;
%[deg_eff,lstr] = legend_info(mx,maxlevel);
sym_color = 'b';

lim = 0;
add2plot = nolim & levels(maxlevel) & info.uniform;
fname = 'level2_uniform_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

lim = 1;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level2_uniform_withlim_1_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

lim = 2;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level2_uniform_withlim_2_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

lim = 3;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform; 
fname = 'level2_uniform_withlim_3_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level2_uniform_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ---------------------
lim = 0;
add2plot = nolim & levels(maxlevel) & info.adapt;
fname = 'level2_adapt_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,20,add2plot);
nf = nf+1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.adapt;
fname = 'level2_adapt_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,8,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=3 ----------------------------------
maxlevel=3;
sym_color = 'g';

lim = 0;
add2plot = nolim & levels(3) & info.uniform;
fname = 'level3_uniform_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

lim = 1;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level3_uniform_withlim_1_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

lim = 2;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level3_uniform_withlim_2_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

lim = 3;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level3_uniform_withlim_3_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level3_uniform_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ---------------------
lim = 0;
add2plot = nolim & levels(maxlevel) & info.adapt;
fname = 'level3_adapt_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,20,add2plot);
nf = nf+1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.adapt;
fname = 'level3_adapt_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,8,add2plot);
nf = nf+1;

% ---------------------------- MAXLEVEL=4 ---------------------------------
maxlevel = 4;
sym_color = 'm';

lim = 0;
add2plot = nolim  & levels(maxlevel) & info.uniform;
fname = 'level4_uniform_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf + 1;

lim = 1;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level4_uniform_withlim_1_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf + 1;

lim = 2;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level4_uniform_withlim_2_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf + 1;

lim = 3;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level4_uniform_withlim_3_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf + 1;

lim = 4;
add2plot = withlim(lim) & levels(maxlevel) & info.uniform;
fname = 'level4_uniform_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf + 1;

% ------------- ADAPT ---------------------
lim = 0;
add2plot = nolim  & levels(maxlevel) & info.adapt;
fname = 'level4_adapt_nolim_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,20,add2plot);
nf = nf + 1;

lim = 4;
add2plot = withlim(lim)  & levels(maxlevel) & info.adapt;
fname = 'level4_adapt_withlim_4_t2.5.dat';
lstr = legend_str(mx,maxlevel,lim,info);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,8,add2plot);
nf = nf + 1;

% -------------------------------------------------------------------------

if (nf ~= nfiles_max+1)
    error('Maxinum number of files exceeded');
end


plot([0 1],[100 100],'k-','linewidth',3);
hold on;
plot([0 1],[0 0],'k-','linewidth',3);

p = zeros(nfiles_max,1);
for i = 1:nfiles_max
    if s(i).add2plot
        fn = [s(i).dir,'/',s(i).file];
        fprintf('Plotting results for ''%s''\n',fn);
        d = load(fn);
        if (~plot_last)
            d(end,:) = nan;
        end
        p(i) = ...
            plot(d(:,1),d(:,2),'color',s(i).color,'marker',...
                s(i).sym,'markersize',s(i).size);
            set(p(i),'linestyle',s(i).line_style,'linewidth',1);
        hold on;    
    end
end
set(gca,'xtick',[0:0.1:(1+1e-8)]);
axis([0 1 -5 180]);
grid on

m = find(p ~= 0);
lh = legend(p(m),{s(m).lstr},'location','southwest');
set(lh,'fontsize',12);

title('Filament diagnostic','fontsize',18);
xlabel('Filament threshold','fontsize',16);
set(gca,'fontsize',14);

shg

end

function s = add_to_plot(dir,fname,lstr,n,c,sty,lstyle,ms,add2plot)

s = struct('file',[],'lstr',[],'n',[],'color',[],'sym',[],...
    'size',[],'add2plot',[]);

s.dir = dir;
s.file = fname;
s.lstr = lstr;
s.n = n;
s.color = c;
s.sym = sty;
s.size = ms;
s.add2plot = add2plot;
s.line_style = lstyle;

end

function lstr = legend_str(mx,maxlevel,lim,info)

global show_degs

if (info.adapt)
    astr = 'A';
end
if (info.uniform)
    astr = 'U';
end
limstr = sprintf('%d',lim);

if (show_degs)
    deg_eff = 90/(mx*2^maxlevel);
    lstr = sprintf('%.2f^o (%s,%d)',deg_eff,astr,lim);
else
    deg_eff = (mx*2^maxlevel);
    lstr = sprintf('%d (%s,%d)',deg_eff,astr,limstr);
end

end
