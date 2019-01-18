%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 12;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';
  
ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
% 
        s_fc(:,:,isubj,m,1,ifoi,iblock) = powcorr;
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
        s_fc(:,:,isubj,m,2,ifoi,iblock) = powcorr;

      end
    end
  end
end

s_fc = nanmean(s_fc(:,:,SUBJLIST,:,:,:,:),7);


%%
for ifoi = 1 : 13

[h,~,~,s]=ttest(s_fc(:,:,:,2,1,ifoi),s_fc(:,:,:,1,1,ifoi),'dim',3);
n_p_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./8010;
n_n_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./8010;

[h,~,~,s]=ttest(s_fc(:,:,:,2,2,ifoi),s_fc(:,:,:,1,2,ifoi),'dim',3);
n_p_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./8010;
n_n_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./8010;

[h,~,~,s]=ttest(s_fc(:,:,:,3,1,ifoi),s_fc(:,:,:,1,1,ifoi),'dim',3);
n_p_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./8010;
n_n_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./8010;

[h,~,~,s]=ttest(s_fc(:,:,:,3,2,ifoi),s_fc(:,:,:,1,2,ifoi),'dim',3);
n_p_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./8010;
n_n_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./8010;

end


%%
clear tperm_cnt1_n tperm_cnt1_p tperm_cnt2_n tperm_cnt2_p
clear tperm_res1_n tperm_res1_p tperm_res2_n tperm_res2_p
clear perm_n_p_atx perm_n_n_atx perm_n_n_dpz perm_n_p_dpz

nperm = 20000;
v = 1;
par.subs = 100;
par.allperms = nperm/par.subs;
alpha = 0.05;
alp = 0.05;

for iperm = 1 : par.allperms
  
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')

  perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_p;
  perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt1_p;
    
  perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_n;
  perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt1_n;
      
  perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_p;
  perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt2_p;

  perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_n;
  perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt2_n;

end


%%

for ifoi = 1 : 13
  
  % COMPUTE EFFECT OF DRUGS ON TASK
  tp_cnt1_n(ifoi) = 1-[sum(n_n_atx(ifoi,2)>perm_n_n_atx(:,ifoi,2))]/nperm;
  tp_cnt1_p(ifoi) = 1-[sum(n_p_atx(ifoi,2)>perm_n_p_atx(:,ifoi,2))]/nperm;
% 
  tp_cnt2_n(ifoi) = 1-[sum(n_n_dpz(ifoi,2)>perm_n_n_dpz(:,ifoi,2))]/nperm;
  tp_cnt2_p(ifoi) = 1-[sum(n_p_dpz(ifoi,2)>perm_n_p_dpz(:,ifoi,2))]/nperm;
  
  % COMPUTE EFFECT OF DRUGS ON REST
  tp_res1_n(ifoi) = 1-[sum(n_n_atx(ifoi,1)>perm_n_n_atx(:,ifoi,1))]/nperm;
  tp_res1_p(ifoi) = 1-[sum(n_p_atx(ifoi,1)>perm_n_p_atx(:,ifoi,1))]/nperm;

  tp_res2_n(ifoi) = 1-[sum(n_n_dpz(ifoi,1)>perm_n_n_dpz(:,ifoi,1))]/nperm;
  tp_res2_p(ifoi) = 1-[sum(n_p_dpz(ifoi,1)>perm_n_p_dpz(:,ifoi,1))]/nperm;
  
  
%   % COMPUTE CHANGE OF CHANGE
  p_d1_p(ifoi) = 1-sum(abs(n_p_atx(ifoi,1)-n_p_atx(ifoi,2))>abs((perm_n_p_atx(:,ifoi,1)-perm_n_p_atx(:,ifoi,2)./4005)))/nperm;
  p_d1_n(ifoi) = 1-sum(abs(n_n_atx(ifoi,1)-n_n_atx(ifoi,2))>abs((perm_n_n_atx(:,ifoi,1)-perm_n_n_atx(:,ifoi,2)./4005)))/nperm;
% 
  p_d2_p(ifoi) = 1-sum(abs(n_p_dpz(ifoi,1)-n_p_dpz(ifoi,2))>abs((perm_n_p_dpz(:,ifoi,1)-perm_n_p_dpz(:,ifoi,2)./4005)))/nperm;
  p_d2_n(ifoi) = 1-sum(abs(n_n_dpz(ifoi,1)-n_n_dpz(ifoi,2))>abs((perm_n_n_dpz(:,ifoi,1)-perm_n_n_dpz(:,ifoi,2)./4005)))/nperm;
% 
  d1_p(ifoi) = n_p_atx(ifoi,1)-n_p_atx(ifoi,2);
  d1_n(ifoi) = n_n_atx(ifoi,1)-n_n_atx(ifoi,2);
  
  d2_p(ifoi) = n_p_dpz(ifoi,1)-n_p_dpz(ifoi,2);
  d2_n(ifoi) = n_n_dpz(ifoi,1)-n_n_dpz(ifoi,2);
  
end


%% PLOT P-VALUES
figure;

subplot(3,2,1); hold on
plot(-log10(tp_res1_n),'b-','linewidth',3)
plot(-log10(tp_res1_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Rest')

subplot(3,2,2); hold on
plot(-log10(tp_res2_n),'b-','linewidth',3)
plot(-log10(tp_res2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

subplot(3,2,3); hold on
plot(-log10(tp_cnt1_n),'b-','linewidth',3)
plot(-log10(tp_cnt1_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Task')

subplot(3,2,4); hold on
plot(-log10(tp_cnt2_n),'b-','linewidth',3)
plot(-log10(tp_cnt2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

subplot(3,2,5); hold on
plot(-log10(p_d1_n),'b-','linewidth',3)
plot(-log10(p_d1_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

subplot(3,2,6); hold on
plot(-log10(p_d2_n),'b-','linewidth',3)
plot(-log10(p_d2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp_pval.pdf'));


%% PLOT CHANGES

figure;

subplot(3,2,1); hold on
plot(n_p_atx(:,1),'r-','linewidth',3)
plot(n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',[0 10 20 30])
ylabel('Altered corr. [%]')
title('Rest')
axis([1 14 -0.05 0.3])
plot(find(tp_res1_p<0.025),n_p_atx(find(tp_res1_p<0.025),1),'k.','markersize',30)
plot(find(tp_res1_n<0.025),n_n_atx(find(tp_res1_n<0.025),1),'k.','markersize',30)

subplot(3,2,2); hold on
plot(n_p_dpz(:,1),'r-','linewidth',3)
plot(n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',[0 10 20 30])
axis([1 14 -0.05 0.3])
plot(find(tp_res2_p<0.025),n_p_dpz(find(tp_res2_p<0.025),1),'k.','markersize',30)
plot(find(tp_res2_n<0.025),n_n_dpz(find(tp_res2_n<0.025),1),'k.','markersize',30)

subplot(3,2,3); hold on
plot(n_p_atx(:,2),'r-','linewidth',3)
plot(n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',[0 10 20 30])
ylabel('Altered corr. [%]')
title('Task')
axis([1 14 -0.05 0.3])
plot(find(tp_cnt1_p<0.025),n_p_atx(find(tp_cnt1_p<0.025),2),'k.','markersize',30)
plot(find(tp_cnt1_n<0.025),n_p_atx(find(tp_cnt1_n<0.025),2),'k.','markersize',30)

subplot(3,2,4); hold on
plot(n_p_dpz(:,2),'r-','linewidth',3)
plot(n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',[0 10 20 30])
axis([1 14 -0.05 0.3])
plot(find(tp_cnt2_p<0.025),n_p_dpz(find(tp_cnt2_p<0.025),2),'k.','markersize',30)
plot(find(tp_cnt2_n<0.025),n_n_dpz(find(tp_cnt2_n<0.025),2),'k.','markersize',30)

subplot(3,2,5); hold on
plot(d1_n,'b-','linewidth',3)
plot(d1_p,'r-','linewidth',3)
axis([0 14 -0.3 0.3])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(p_d1_p<0.025),d1_p(find(p_d1_p<0.025)),'k.','markersize',30)
plot(find(p_d1_n<0.025),d1_n(find(p_d1_n<0.025)),'k.','markersize',30)

subplot(3,2,6); hold on
plot(d2_n,'b-','linewidth',3)
plot(d2_p,'r-','linewidth',3)
axis([0 14 -0.3 0.3])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); 
plot(find(p_d2_p<0.025),d2_p(find(p_d2_p<0.025)),'k.','markersize',30)
plot(find(p_d2_n<0.025),d2_n(find(p_d2_n<0.025)),'k.','markersize',30)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp.pdf'));



%% MULTIPLE COMPARISONS CORRECTION (HAWELLEK ET AL., xxxx)

for ifreq = 1 : 13
  [~,idx_D_cnt1_n(:,ifreq)] = sort(perm_n_n_atx(:,ifreq,2));
  [~,idx_D_cnt1_p(:,ifreq)] = sort(perm_n_p_atx(:,ifreq,2));
  [~,idx_D_cnt2_n(:,ifreq)] = sort(perm_n_n_dpz(:,ifreq,2));
  [~,idx_D_cnt2_p(:,ifreq)] = sort(perm_n_p_dpz(:,ifreq,2));
  [~,idx_D_res1_n(:,ifreq)] = sort(perm_n_n_atx(:,ifreq,1));
  [~,idx_D_res1_p(:,ifreq)] = sort(perm_n_p_atx(:,ifreq,1));
  [~,idx_D_res2_n(:,ifreq)] = sort(perm_n_n_dpz(:,ifreq,1));
  [~,idx_D_res2_p(:,ifreq)] = sort(perm_n_p_dpz(:,ifreq,1));
end

idx_R_cnt1_p = max(idx_D_cnt1_p,[],2);
idx_R_cnt1_n = max(idx_D_cnt1_n,[],2);
idx_R_cnt2_p = max(idx_D_cnt2_p,[],2);
idx_R_cnt2_n = max(idx_D_cnt2_n,[],2);
idx_R_res1_n = max(idx_D_res1_n,[],2);
idx_R_res1_p = max(idx_D_res1_p,[],2);
idx_R_res2_n = max(idx_D_res2_n,[],2);
idx_R_res2_p = max(idx_D_res2_p,[],2);

for ifreq = 1 : 13
  for irank = 1 : length(idx_R_cnt1_p)
    irank

    tmp = perm_n_p_atx(find(idx_D_res1_p(:,ifreq) == idx_R_res1_p(irank)),ifreq,1);
    tperm_res1_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_res1_p(:,ifreq) == irank)),ifreq,1);
    tmp = perm_n_n_atx(find(idx_D_res1_n(:,ifreq) == idx_R_res1_n(irank)),ifreq,1);
    tperm_res1_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_res1_n(:,ifreq) == irank)),ifreq,1);
  
    tmp = perm_n_p_dpz(find(idx_D_res2_p(:,ifreq) == idx_R_res2_p(irank)),ifreq,1);
    tperm_res2_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_res2_p(:,ifreq) == irank)),ifreq,1);
    tmp = perm_n_n_dpz(find(idx_D_res2_n(:,ifreq) == idx_R_res2_n(irank)),ifreq,1);
    tperm_res2_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_res2_n(:,ifreq) == irank)),ifreq,1);

    tmp = perm_n_p_atx(find(idx_D_cnt1_p(:,ifreq) == idx_R_cnt1_p(irank)),ifreq,2);
    tperm_cnt1_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_cnt1_p(:,ifreq) == irank)),ifreq,2);
    tmp = perm_n_n_atx(find(idx_D_cnt1_n(:,ifreq) == idx_R_cnt1_n(irank)),ifreq,2);
    tperm_cnt1_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_cnt1_n(:,ifreq) == irank)),ifreq,2);
     
    tmp = perm_n_p_dpz(find(idx_D_cnt2_p(:,ifreq) == idx_R_cnt2_p(irank)),ifreq,2);
    tperm_cnt2_p_corr(irank,ifreq)=perm_n_p_dpz(unique(find(idx_D_cnt2_p(:,ifreq) == irank)),ifreq,2);   
    tmp = perm_n_n_dpz(find(idx_D_cnt2_n(:,ifreq) == idx_R_cnt2_n(irank)),ifreq,2);
    tperm_cnt2_n_corr(irank,ifreq)=perm_n_n_dpz(unique(find(idx_D_cnt2_n(:,ifreq) == irank)),ifreq,2);
  

  end
end

for ifreq = 1 : 13

  p_res1_p(ifreq) = 1-sum(n_p_atx(ifreq,1)>tperm_res1_p_corr(:,ifreq))/nperm;
  p_res1_n(ifreq) = 1-sum(n_n_atx(ifreq,1)>tperm_res1_n_corr(:,ifreq))/nperm;
  p_cnt1_p(ifreq) = 1-sum(n_p_atx(ifreq,2)>tperm_cnt1_p_corr(:,ifreq))/nperm;
  p_cnt1_n(ifreq) = 1-sum(n_n_atx(ifreq,2)>tperm_cnt1_n_corr(:,ifreq))/nperm;
  p_res2_p(ifreq) = 1-sum(n_p_dpz(ifreq,1)>tperm_res2_p_corr(:,ifreq))/nperm;
  p_res2_n(ifreq) = 1-sum(n_n_dpz(ifreq,1)>tperm_res2_n_corr(:,ifreq))/nperm;
  p_cnt2_p(ifreq) = 1-sum(n_p_dpz(ifreq,2)>tperm_cnt2_p_corr(:,ifreq))/nperm;
  p_cnt2_n(ifreq) = 1-sum(n_n_dpz(ifreq,2)>tperm_cnt2_n_corr(:,ifreq))/nperm;

end


%% BEHAV

m = squeeze(s_fc(:,:,:,2,2,6))-squeeze(s_fc(:,:,:,1,2,6));

m = squeeze(nanmean(nanmean(m),2))

s = 1:28;

% for isubj = 1 : 28
%   
%   ss = s(s~=isubj);
%   
%   p_fc(isubj) = sum(sum(triu(1-(sum((s_fc(:,:,isubj,2,2,6)-s_fc(:,:,isubj,1,2,6)) > m,3)./27)<0.05,1)))./4005;
%   n_fc(isubj) = sum(sum(triu(1-(sum((s_fc(:,:,isubj,2,2,6)-s_fc(:,:,isubj,1,2,6)) < m,3)./27)<0.05,1)))./4005;
%   
% end



