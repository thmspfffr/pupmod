%% LOAD AND PLOT FRACTION OF NODES WITH ALTERED VARIANCE
% ------------------------------------

clear

outdir = '~/pupmod/proc/conn/';

% Version 23: 400 vertices, cortical surface
v =3;

% Include all 28 paricipants
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

ord   = pconn_randomization;
for m = 1 : 3
  for isubj = SUBJLIST
    isubj
        
    im = find(ord(isubj,:)==m);
    
    load(sprintf([outdir 'pupmod_src_variance_s%d_m%d_v%d.mat'],isubj,m,v));
    
    var_all_rest(:,isubj,m,:,:) = variance; 
    
    clear variance
    
    load(sprintf([outdir 'pupmod_task_src_variance_s%d_m%d_v%d.mat'],isubj,m,v));
    
    var_all_task(:,isubj,m,:,:) = variance; 
    
    clear variance

  end
end

% average across recording blocks
var_all_rest = squeeze(nanmean(var_all_rest(:,SUBJLIST,:,:,:),4));
var_all_task = squeeze(nanmean(var_all_task(:,SUBJLIST,:,:,:),4));

%% PERMUTATION TEST
% Run permutation test, without correction for multiple comparisons
% (across frequencies). N(perm.)=10.000

for ifoi = 1 : 17
  
  % index 2 = atx, index 1 = pbo
  [h,~,~,s]=ttest(var_all_rest(:,:,2,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_atx_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  % index 3 = dpz, index 1 = pbo
  [h,~,~,s]=ttest(var_all_rest(:,:,3,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_dpz_pos_rest(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_dpz_neg_rest(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  % index 2 = atx, index 1 = pbo
  [h,~,~,s]=ttest(var_all_task(:,:,2,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
  n_atx_pos_task(ifoi) = 100*sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_task(ifoi) = 100*sum((h>0)&(s.tstat<0))./ length(h);
  
  % index 3 = dpz, index 1 = pbo
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
  
  for ifoi = 1 : 17
    
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

load(sprintf([outdir 'pupmod_src_variance_permutation_v%d.mat'],v));

% ATOMOXETINE
% ------------------
idx_R_res1_p      = max(abs(perm.n_atx_pos_rest),[],2);
idx_R_res1_n      = max(abs(perm.n_atx_neg_rest),[],2);
idx_R_cnt1_p      = max(abs(perm.n_atx_pos_task),[],2);
idx_R_cnt1_n      = max(abs(perm.n_atx_neg_task),[],2);

% DONEPEZIL
% ------------------
idx_R_res2_p      = max(abs(perm.n_dpz_pos_rest),[],2);
idx_R_res2_n      = max(abs(perm.n_dpz_neg_rest),[],2);
idx_R_cnt2_p      = max(abs(perm.n_dpz_pos_task),[],2);
idx_R_cnt2_n      = max(abs(perm.n_dpz_neg_task),[],2);


for ifoi = 1:17
  % ATOMOXETINE
  % ------------------
  p_atx_pos_rest(ifoi) = 1-sum(abs(n_atx_pos_rest(ifoi))>abs(idx_R_res1_p))/nperm;
  p_atx_neg_rest(ifoi) = 1-sum(abs(n_atx_neg_rest(ifoi))>abs(idx_R_res1_n))/nperm;
  p_atx_pos_task(ifoi) = 1-sum(abs(n_atx_pos_task(ifoi))>abs(idx_R_cnt1_p))/nperm;
  p_atx_neg_task(ifoi) = 1-sum(abs(n_atx_neg_task(ifoi))>abs(idx_R_cnt1_n))/nperm;

  p_dpz_pos_rest(ifoi) = 1-sum(abs(n_dpz_pos_rest(ifoi))>abs(idx_R_res2_p))/nperm;
  p_dpz_neg_rest(ifoi) = 1-sum(abs(n_dpz_neg_rest(ifoi))>abs(idx_R_res2_n))/nperm;
  p_dpz_pos_task(ifoi) = 1-sum(abs(n_dpz_pos_task(ifoi))>abs(idx_R_cnt2_p))/nperm;
  p_dpz_neg_task(ifoi) = 1-sum(abs(n_dpz_neg_task(ifoi))>abs(idx_R_cnt2_n))/nperm;

end


%%
foi_range       = 2.^[2:.25:6];

markersize = 2;
figure; set(gcf,'color','w')
% -----------
% ATOMOXETINE DURING REST
% -----------
subplot(5,4,1); hold on

plot(n_atx_pos_rest,'k-','linewidth',0.75)
plot(n_atx_neg_rest,':','linewidth',0.75,'color',[.44 .44 .44])

plot(find(p_atx_pos_rest<0.05),n_atx_pos_rest(find(p_atx_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_atx_neg_rest<0.05),n_atx_neg_rest(find(p_atx_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 17 -2 25]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17 ],'xticklabel',[4 8 16 32 64])
set(gca,'tickdir','out','ytick',[0 25],'yticklabel',num2cell([0 25]))

ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
tp_editplots
% -----------
% ATOMOXETINE DURING TASK
% -----------
subplot(5,4,5); hold on
plot(n_atx_pos_task,'k-','linewidth',0.75)
plot(n_atx_neg_task,':','linewidth',0.75,'color',[.44 .44 .44])

plot(find(p_atx_pos_task<0.05),n_atx_pos_task(find(p_atx_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_atx_neg_task<0.05),n_atx_neg_task(find(p_atx_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 17 -2 25]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
set(gca,'tickdir','out','ytick',[0 25],'yticklabel',num2cell([0 25]))

xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
tp_editplots
% -----------
% DONEPEZIL DURING REST
% -----------
subplot(5,4,2); hold on
plot(n_dpz_pos_rest,'k-','linewidth',0.75)
plot(n_dpz_neg_rest,':','linewidth',0.75,'color',[.44 .44 .44])

plot(find(p_dpz_pos_rest<0.05),n_dpz_pos_rest(find(p_dpz_pos_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_dpz_neg_rest<0.05),n_dpz_neg_rest(find(p_dpz_neg_rest<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 17 -2 25]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
set(gca,'tickdir','out','ytick',[0 25],'yticklabel',num2cell([0 25]))

tp_editplots

% -----------
% DONEPEZIL DURING TASK
% -----------
subplot(5,4,6); hold on
plot(n_dpz_pos_task,'k-','linewidth',0.75)
plot(n_dpz_neg_task,':','linewidth',0.75,'color',[.44 .44 .44])
plot(find(p_dpz_pos_task<0.05),n_dpz_pos_task(find(p_dpz_pos_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(p_dpz_neg_task<0.05),n_dpz_neg_task(find(p_dpz_neg_task<0.05)),'ko','markersize',markersize,'markerfacecolor','k')

axis([0 17 -2 25]);
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
set(gca,'tickdir','out','ytick',[0 25],'yticklabel',num2cell([0 25]))
xlabel('Carrier frequency [Hz]'); 
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_variance_v%d.eps',v))



