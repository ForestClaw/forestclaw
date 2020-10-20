function plot_errors(e,idx,compare)

% l = size(e,1);

if (nargin < 3)    
    compare = false;
    if (nargin < 2)
        n = size(e,1);
        idx = 1:n;
    end
end

lh = zeros(3,1);
lsq = zeros(3,1);
lstr_basic = {'1-norm','2-norm','inf-norm'};
lstr = lstr_basic;
if (~compare)
    close all
    ladd = 0;
else
    lh_compare = findobj('Tag','compare');
    lsq = findobj('Tag','lsqfit');
    if (~isempty(lh_compare))        
        lh = lh_compare;
        o = findobj('Tag','legend');
        lstr = o.String;
        ladd = 3;
    else
        ladd = 0;
    end
    if (~isempty(lsq))
        delete(lsq);
    end
end

Nvec = e(:,1);
n = length(Nvec);

i1 = min(idx);
if (i1 >= n)
    error('min(idx) >= length(e)');
end

if length(idx) == 1
    idx = i1:length(e);
end
i2 = min([max(idx),n]);
idx = i1:i2;

c = {'r','b','g'};
if (nargin == 1)
    idx = 1:n;
end
for i = 1:3
    p = polyfit(log(Nvec(idx)),log(e(idx,i+1)),1);
    if (~compare)
        lsq(i) = loglog(Nvec,exp(polyval(p,log(Nvec))),'k','linewidth',1); 
        set(lsq(i),'Tag','lsqfit');
        hold on;
    end
    lh(i+ladd) = loglog(Nvec,e(:,i+1),c{i},'linewidth',1);
    lstr{i+ladd} = sprintf('%10s (rate = %6.4f)',lstr_basic{i},-p(1));
    hold on;
end

if (~compare)
    mh = loglog(Nvec,e(:,2:4),'ko','markersize',8);
    mh_idx = loglog(Nvec(idx),e(idx,2:4),'k.','markersize',28);
else
    n = length(lh);
    set(lh(end-2:end),'linestyle','--');
    mh = loglog(Nvec,e(:,2:4),'kp','markersize',12);
    % set(mh,'markerfacecolor','r');
    mh_idx = loglog(Nvec(idx),e(idx,2:4),'kp','markersize',14);
    set(mh_idx,'markerfacecolor','k');
end    
set(lh,'Tag','compare');

legend(lh,lstr,'location','southwest');

xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);
title('Error','fontsize',18);

set(gca,'xtick',Nvec);

p0 = log2(Nvec(1));
p1 = log2(Nvec(end));
xlim([2^(p0-0.5), 2^(p1+0.5)]);

set(gca,'fontsize',16);

m2 = size(e,2);
idx = 2:m2;
        
conv_rates = log(e(1:end-1,idx)./e(2:end,idx))/log(2);
cr_data = [Nvec(2:end) Nvec(1:end-1) conv_rates]';

% Create variable lenght print string
cr_cols = size(conv_rates,2);
cr_str = kron(ones(1,cr_cols),double('%8.4f '));
print_str = sprintf('%s %s\\n','%5d/%5d',cr_str); 

% Draw line
cr_len = 12 + 9*cr_cols;
cr_line = double('-')*ones(1,cr_len);

% Print table
fprintf('\n');
fprintf('          Convergence rates\n');
fprintf('%s\n',cr_line);
fprintf(print_str,cr_data);
fprintf('%s\n',cr_line);
    
end