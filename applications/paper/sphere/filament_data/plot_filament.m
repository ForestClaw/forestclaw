function plot_filament(dir)

close all

% Limiters
% 0 - no limiter
% 1 - Minmod
% 2 - Superbee
% 3 - Van Leer
% 4 - Montonized Centered

nolim = false;
withlim = [false, false, false, true];  % Use limiter ID
levels = [true, true, true, true];
uniform = false;
adapt = true;

lsu = ':';  % line style for uniform
lsa = '-'; % line style for adaptive

mx = 30;
nf = 1;    % file counter (incremented each time a file is added)

nfiles_max = 30;   % Maximum number of files
% ------------------------- Initial conditions ----------------------------
fname = 'level2_uniform_nolim_t0.dat';
lstr = sprintf('t=0.0');
s = add_to_plot(dir,fname,lstr,nf,'r','.',lsu,30,false);
nf = nf+1;

% -------------------------------- Test -----------------------------------
add2plot = false;
fname = 'filament.dat';
lstr = sprintf('TEST');
s = add_to_plot(dir,fname,lstr,nf,'k','.','-',30,add2plot);
nf = nf+1;

% ----------------------------- MAXLEVEL=1 --------------------------------
maxlevel=1;
deg_eff = 90/(mx*2^maxlevel);
sym_color = 'r';

add2plot = nolim & levels(1) & uniform;
fname = 'level1_uniform_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (U,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

add2plot = withlim(1) & levels(1) & uniform;
fname = 'level1_uniform_withlim_1_t2.5.dat';
lstr = sprintf('%.2f^o (U,1)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(2) & levels(1) & uniform;
fname = 'level1_uniform_withlim_2_t2.5.dat';
lstr = sprintf('%.2f^o (U,2)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(3) & levels(1)  & uniform; 
fname = 'level1_uniform_withlim_3_t2.5.dat';
lstr = sprintf('%.2f^o (U,3)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(1) & uniform;
fname = 'level1_uniform_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (U,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ---------------------
add2plot = nolim & levels(1) & adapt;
fname = 'level1_adapt_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (A,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,30,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(1) & adapt;
fname = 'level1_adapt_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (A,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,10,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=2 ----------------------------------
maxlevel=2;
deg_eff = 90/(mx*2^maxlevel);
sym_color = 'b';

add2plot = nolim & levels(2) & uniform;
fname = 'level2_uniform_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (U,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

add2plot = withlim(1) & levels(2) & uniform;
fname = 'level2_uniform_withlim_1_t2.5.dat';
lstr = sprintf('%.2f^o (U,1)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(2) & levels(2) & uniform;
fname = 'level2_uniform_withlim_2_t2.5.dat';
lstr = sprintf('%.2f^o (U,2)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(3) & levels(2) & uniform; 
fname = 'level2_uniform_withlim_3_t2.5.dat';
lstr = sprintf('%.2f^o (U,3)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(2) & uniform;
fname = 'level2_uniform_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (U,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ---------------------
add2plot = nolim & levels(2) & adapt;
fname = 'level2_adapt_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (A,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,30,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(2) & adapt;
fname = 'level2_adapt_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (A,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,10,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=3 ----------------------------------
maxlevel=3;
deg_eff = 90/(mx*2^maxlevel);
sym_color = 'g';

add2plot = nolim & levels(3) & uniform;
fname = 'level3_uniform_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (U,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf+1;

add2plot = withlim(1) & levels(3) & uniform;
fname = 'level3_uniform_withlim_1_t2.5.dat';
lstr = sprintf('%.2f^o (U,1)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(2) & levels(3) & uniform;
fname = 'level3_uniform_withlim_2_t2.5.dat';
lstr = sprintf('%.2f^o (U,2)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(3) & levels(3) & uniform;
fname = 'level3_uniform_withlim_3_t2.5.dat';
lstr = sprintf('%.2f^o (U,3)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(3) & uniform;
fname = 'level3_uniform_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (U,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf+1;

% ------------- ADAPT ---------------------
add2plot = nolim & levels(3) & adapt;
fname = 'level3_adapt_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (A,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,30,add2plot);
nf = nf+1;

add2plot = nolim & levels(3) & adapt;
fname = 'level3_adapt_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (A,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,10,add2plot);
nf = nf+1;

% ---------------------------- MAXLEVEL=4 ---------------------------------
maxlevel=4;
deg_eff = 90/(mx*2^maxlevel);
sym_color = 'm';

add2plot = nolim  & levels(4) & uniform;
fname = 'level4_uniform_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (U,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsu,30,add2plot);
nf = nf + 1;

add2plot = withlim(1) & levels(4) & uniform;
fname = 'level4_uniform_withlim_1_t2.5.dat';
lstr = sprintf('%.2f^o (U,1)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsu,10,add2plot);
nf = nf + 1;

add2plot = withlim(2) & levels(4) & uniform;
fname = 'level4_uniform_withlim_2_t2.5.dat';
lstr = sprintf('%.2f^o (U,2)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'d',lsu,10,add2plot);
nf = nf + 1;

add2plot = withlim(3) & levels(4) & uniform;
fname = 'level4_uniform_withlim_3_t2.5.dat';
lstr = sprintf('%.2f^o (U,3)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'*',lsu,10,add2plot);
nf = nf + 1;

add2plot = withlim(4) & levels(4) & uniform;
fname = 'level4_uniform_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (U,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsu,10,add2plot);
nf = nf + 1;

% ------------- ADAPT ---------------------
add2plot = nolim  & levels(4) & adapt;
fname = 'level4_adapt_nolim_t2.5.dat';
lstr = sprintf('%.2f^o (A,0)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,30,add2plot);
nf = nf + 1;

add2plot = withlim(4)  & levels(4) & adapt;
fname = 'level4_adapt_withlim_4_t2.5.dat';
lstr = sprintf('%.2f^o (A,4)',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,10,add2plot);
nf = nf + 1;

% -------------------------------------------------------------------------

if (nf > nfiles_max+1)
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
        p(i) = ...
            plot(d(:,1),d(:,2),'color',s(i).color,'marker',...
                s(i).sym,'markersize',s(i).size);
            set(p(i),'linestyle',s(i).line_style,'linewidth',1);
        hold on;    
    end
end
set(gca,'xtick',[0:0.1:(1+1e-8)]);
axis([0.1 1 -5 140]);
grid on

m = find(p ~= 0);
lh = legend(p(m),{s(m).lstr},'location','southwest');
set(lh,'fontsize',12);

title('Filament diagnostic (cosine hill)','fontsize',18);
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
