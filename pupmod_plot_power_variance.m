%% LOAD VARIANCE DATA
outdir = '~/pupmod/proc/conn/'
clear var_all_rest
clear var_all_task
v =23;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
ord   = pconn_randomization;
for m = 1 : 3
  for isubj = SUBJLIST
    isubj
        
    im = find(ord(isubj,:)==m);
    
    load(sprintf([outdir 'pupmod_src_power_variance_s%d_m%d_v%d.mat'],isubj,m,v));
    
    var_all_rest(:,isubj,m,:,:) = outp.var; 
    pow_all_rest(:,isubj,m,:,:) = outp.pow; 

    load(sprintf([outdir 'pupmod_src_power_variance_task_s%d_m%d_v%d.mat'],isubj,m,v));
    
    var_all_task(:,isubj,m,:,:) = outp.var; 
    pow_all_task(:,isubj,m,:,:) = outp.pow; 

  end
end

var_all_rest = nanmean(var_all_rest(:,SUBJLIST,:,:,:),5);
var_all_task = nanmean(var_all_task(:,SUBJLIST,:,:,:),5);
pow_all_rest = nanmean(pow_all_rest(:,SUBJLIST,:,:,:),5);
pow_all_task = nanmean(pow_all_task(:,SUBJLIST,:,:,:),5);

%% EMPIRICAL

foi_range = [2 3 4 6 8 11 16 23 32 45 64 91 128];

for ifoi = 1 : 25
  
  [h,~,~,s]=ttest(var_all_rest(:,:,2,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_atx_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_rest(:,:,3,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_dpz_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_dpz_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_task(:,:,2,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
  n_atx_pos_task(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_task(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_task(:,:,3,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
  n_dpz_pos_task(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_dpz_neg_task(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
end

% PERMUTATION
nperm = 10000;
all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);

dat_cnt1 = var_all_task(:,:,[1 2],:);
dat_res1 = var_all_rest(:,:,[1 2],:);
dat_cnt2 = var_all_task(:,:,[1 3],:);
dat_res2 = var_all_rest(:,:,[1 3],:);

for iperm = 1 : nperm
  
  % within subjects permutation test
  fprintf('Perm #%d ...\n',iperm);
  
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    
    permdat_cnt1(:,i,1,:) = dat_cnt1(:,i,idx1(i),:);
    permdat_cnt1(:,i,2,:) = dat_cnt1(:,i,idx2(i),:);
    
    permdat_res1(:,i,1,:) = dat_res1(:,i,idx1(i),:);
    permdat_res1(:,i,2,:) = dat_res1(:,i,idx2(i),:);
    
  end
  
  for i = 1 : length(idx1)
    
    permdat_cnt2(:,i,1,:) = dat_cnt2(:,i,idx1(i),:);
    permdat_cnt2(:,i,2,:) = dat_cnt2(:,i,idx2(i),:);
    
    permdat_res2(:,i,1,:) = dat_res2(:,i,idx1(i),:);
    permdat_res2(:,i,2,:) = dat_res2(:,i,idx2(i),:);
    
  end
  
  for ifoi = 1 : 25
    
    [h,~,~,s]=ttest(permdat_res1(:,:,2,ifoi),permdat_res1(:,:,1,ifoi),'dim',2);
    perm.n_atx_pos_rest(iperm,ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_atx_neg_rest(iperm,ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(permdat_res2(:,:,2,ifoi),permdat_res2(:,:,1,ifoi),'dim',2);
    perm.n_dpz_pos_rest(iperm,ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_dpz_neg_rest(iperm,ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(permdat_cnt1(:,:,2,ifoi),permdat_cnt1(:,:,1,ifoi),'dim',2);
    perm.n_atx_pos_task(iperm,ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_atx_neg_task(iperm,ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(permdat_cnt2(:,:,2,ifoi),permdat_cnt2(:,:,1,ifoi),'dim',2);
    perm.n_dpz_pos_task(iperm,ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_dpz_neg_task(iperm,ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
    
    
  end
end

save(sprintf([outdir 'pupmod_src_variance_permutation_v%d.mat'],v),'perm');

for ifoi = 1 : 25
  
  p_atx_pos_rest(ifoi) = 1-(sum(abs(n_atx_pos_rest(:,ifoi))>abs(perm.n_atx_pos_rest(:,ifoi)))./nperm);
  p_atx_neg_rest(ifoi) = 1-(sum(abs(n_atx_neg_rest(:,ifoi))>abs(perm.n_atx_neg_rest(:,ifoi)))./nperm);
  p_atx_pos_task(ifoi) = 1-(sum(abs(n_atx_pos_task(:,ifoi))>abs(perm.n_atx_pos_task(:,ifoi)))./nperm);
  p_atx_neg_task(ifoi) = 1-(sum(abs(n_atx_neg_task(:,ifoi))>abs(perm.n_atx_neg_task(:,ifoi)))./nperm);
  
  p_dpz_pos_rest(ifoi) = 1-(sum(abs(n_dpz_pos_rest(:,ifoi))>abs(perm.n_dpz_pos_rest(:,ifoi)))./nperm);
  p_dpz_neg_rest(ifoi) = 1-(sum(abs(n_dpz_neg_rest(:,ifoi))>abs(perm.n_dpz_neg_rest(:,ifoi)))./nperm);
  p_dpz_pos_task(ifoi) = 1-(sum(abs(n_dpz_pos_task(:,ifoi))>abs(perm.n_dpz_pos_task(:,ifoi)))./nperm);
  p_dpz_neg_task(ifoi) = 1-(sum(abs(n_dpz_neg_task(:,ifoi))>abs(perm.n_dpz_neg_task(:,ifoi)))./nperm);
  
end

%%

markersize = 5;
figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(n_atx_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_atx_neg_rest,'linewidth',2,'color',[0.2 0.5 1])

plot(find(p_atx_pos_rest<0.05),n_atx_pos_rest(find(p_atx_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_atx_neg_rest<0.05),n_atx_neg_rest(find(p_atx_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 21 -2 30]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Frequency [Hz]');
ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,3); hold on
plot(n_atx_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_atx_neg_task,'linewidth',2,'color',[0.2 0.5 1])

plot(find(p_atx_pos_task<0.05),n_atx_pos_task(find(p_atx_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_atx_neg_task<0.05),n_atx_neg_task(find(p_atx_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

% plot(prctile(n_atx_task_perm,95),'linewidth',1,'color',[1 0.1 0.1],'linestyle',':')
axis([0 21 -2 30]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,2); hold on
plot(n_dpz_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_dpz_neg_rest,'linewidth',2,'color',[0.2 0.5 1])

plot(find(p_dpz_pos_rest<0.05),n_dpz_pos_rest(find(p_dpz_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_dpz_neg_rest<0.05),n_dpz_neg_rest(find(p_dpz_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 21 -2 30]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Donepezil vs. placebo')
tp_editplots

subplot(4,2,4); hold on
plot(n_dpz_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_dpz_neg_task,'linewidth',2,'color',[0.2 0.5 1])
plot(find(p_dpz_pos_task<0.05),n_dpz_pos_task(find(p_dpz_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_dpz_neg_task<0.05),n_dpz_neg_task(find(p_dpz_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

% plot(prctile(n_dpz_task_perm,95),'linewidth',1,'color',[0.1 0.1 1],'linestyle',':')
axis([0 21 -2 30]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 10 20 30],'yticklabel',num2cell([0 10 20 30]))
xlabel('Carrier frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_variance_v%d.eps',v))
%% POWER 


for ifoi = 1 : 25
  
  [h,~,~,s]=ttest(pow_all_rest(:,:,2,ifoi),pow_all_rest(:,:,1,ifoi),'dim',2);
  n_pow_atx_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_pow_atx_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(pow_all_rest(:,:,3,ifoi),pow_all_rest(:,:,1,ifoi),'dim',2);
  n_pow_dpz_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_pow_dpz_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(pow_all_task(:,:,2,ifoi),pow_all_task(:,:,1,ifoi),'dim',2);
  n_pow_atx_pos_task(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_pow_atx_neg_task(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(pow_all_task(:,:,3,ifoi),pow_all_task(:,:,1,ifoi),'dim',2);
  n_pow_dpz_pos_task(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_pow_dpz_neg_task(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
end

%%
markersize = 5;
figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(n_pow_atx_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_pow_atx_neg_rest,'linewidth',2,'color',[0.2 0.5 1])

% plot(find(p_atx_pos_rest<0.05),n_atx_pos_rest(find(p_atx_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
% plot(find(p_atx_neg_rest<0.05),n_atx_neg_rest(find(p_atx_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 21 -2 60]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Frequency [Hz]');
ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,3); hold on
plot(n_pow_atx_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_pow_atx_neg_task,'linewidth',2,'color',[0.2 0.5 1])

% plot(find(p_atx_pos_task<0.05),n_atx_pos_task(find(p_atx_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
% plot(find(p_atx_neg_task<0.05),n_atx_neg_task(find(p_atx_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

% plot(prctile(n_atx_task_perm,95),'linewidth',1,'color',[1 0.1 0.1],'linestyle',':')
axis([0 21 -2 60]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,2); hold on
plot(n_pow_dpz_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_pow_dpz_neg_rest,'linewidth',2,'color',[0.2 0.5 1])

% plot(find(p_dpz_pos_rest<0.05),n_dpz_pos_rest(find(p_dpz_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
% plot(find(p_dpz_neg_rest<0.05),n_dpz_neg_rest(find(p_dpz_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 21 -2 60]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Donepezil vs. placebo')
tp_editplots

subplot(4,2,4); hold on
plot(n_pow_dpz_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_pow_dpz_neg_task,'linewidth',2,'color',[0.2 0.5 1])
% plot(find(p_dpz_pos_task<0.05),n_dpz_pos_task(find(p_dpz_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
% plot(find(p_dpz_neg_task<0.05),n_dpz_neg_task(find(p_dpz_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

% plot(prctile(n_dpz_task_perm,95),'linewidth',1,'color',[0.1 0.1 1],'linestyle',':')
axis([0 21 -2 60]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 10 20 30],'yticklabel',num2cell([0 10 20 30]))
xlabel('Carrier frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_power_v%d.eps',v))






