function plot_filament_short(dir)

close all

% Limiters
% 0 - no limiter
% 1 - Minmod
% 2 - Superbee
% 3 - Van Leer
% 4 - Montonized Centered

nolim   =  false;
withlim = [false, false, false, true];     % Use limiter ID
levels = [true, true, true, true, true];      % Levels 1,2,3,4
adapt = true;

lsu = ':';  % line style for uniform
lsa = '-'; % line style for adaptive

mx = 32;
nf = 1;    % file counter (incremented each time a file is added)

nfiles_max = 12;   % Maximum number of files
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
deg_eff = (mx*2^maxlevel);
sym_color = 'r';

add2plot = nolim & levels(1) & adapt;
fname = 'level1_adapt_nolim_t2.5.dat';
% lstr = sprintf('%.2f^o (A,0)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,20,add2plot);
nf = nf+1;

% Adaptive runs only done for limiter=4
add2plot = withlim(4) & levels(1) & adapt;
fname = 'level1_adapt_withlim_4_t2.5.dat';
% lstr = sprintf('%.2f^o (A,4)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'.',lsa,20,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=2 ----------------------------------
maxlevel=2;
deg_eff = 90/(mx*2^maxlevel);
deg_eff = (mx*2^maxlevel);
sym_color = 'b';

add2plot = nolim & levels(2) & adapt;
fname = 'level2_adapt_nolim_t2.5.dat';
% lstr = sprintf('%.2f^o (A,0)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'o',lsa,20,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(2) & adapt;
fname = 'level2_adapt_withlim_4_t2.5.dat';
%lstr = sprintf('%.2f^o (A,4)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'o',lsa,8,add2plot);
nf = nf+1;

% --------------------------- MAXLEVEL=3 ----------------------------------
maxlevel=3;
deg_eff = 90/(mx*2^maxlevel);
deg_eff = (mx*2^maxlevel);
sym_color = 'k';

add2plot = nolim & levels(3) & adapt;
fname = 'level3_adapt_nolim_t2.5.dat';
%lstr = sprintf('%.2f^o (A,0)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'p',lsa,10,add2plot);
nf = nf+1;

add2plot = withlim(4) & levels(3) & adapt;
fname = 'level3_adapt_withlim_4_t2.5.dat';
%lstr = sprintf('%.2f^o (A,4)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'p',lsa,10,add2plot);
nf = nf+1;

% ---------------------------- MAXLEVEL=4 ---------------------------------
maxlevel=4;
deg_eff = (mx*2^maxlevel);
sym_color = 'm';

add2plot = nolim  & levels(4) & adapt;
fname = 'level4_adapt_nolim_t2.5.dat';
%lstr = sprintf('%.2f^o (A,0)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsa,20,add2plot);
nf = nf + 1;

add2plot = withlim(4)  & levels(4) & adapt;
fname = 'level4_adapt_withlim_4_t2.5.dat';
%lstr = sprintf('%.2f^o (A,4)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'s',lsa,8,add2plot);
nf = nf + 1;

% ---------------------------- MAXLEVEL=5 ---------------------------------
maxlevel=5;
deg_eff = (mx*2^maxlevel);
sym_color = 'g';

add2plot = nolim  & levels(5) & adapt;
fname = 'level5_adapt_nolim_t2.5.dat';
%lstr = sprintf('%.2f^o (A,0)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,20,add2plot);
nf = nf + 1;

add2plot = withlim(4)  & levels(5) & adapt;
fname = 'level5_adapt_withlim_4_t2.5.dat';
%lstr = sprintf('%.2f^o (A,4)',deg_eff);
lstr = sprintf('%d',deg_eff);
s(nf) = add_to_plot(dir,fname,lstr,nf,sym_color,'v',lsa,8,add2plot);
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
        d([end],:) = nan;
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

t = 2.5;
tstr = sprintf('Filament diagnostic (t=%.1f)',t);
title(tstr,'fontsize',18);
xlabel('Filament threshold','fontsize',16);
ylabel('A(\tau,2.5)/A(\tau,0) (%)','fontsize',16);
set(gca,'fontsize',14);
grid on

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
