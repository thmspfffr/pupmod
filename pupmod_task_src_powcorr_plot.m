%% pupmod_src_powcorr_plot

clear all

% --------------------------------------------------------
% VERSION 10
% --------------------------------------------------------
v               = 18;
v_postproc      = 6;
v_grid          = 4; % 4 = aal
fsample         = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal';
foi_range       = unique(round(2.^[1:.5:7]));
para.smo        = foi_range./4;
para.segleng    = 1 ./ para.smo;
lpc             = 0;
timevariant     = 1;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/


%%
clear s s1 s2 fc_mean
v = 10;


addpath ~/pconn/matlab/
  
ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        s_cnt(:,:,isubj,m,ifoi,iblock) =  powcorr;
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        s_res(:,:,isubj,m,ifoi,iblock) =  powcorr;
      
      end
    end
  end
end

% save('~/pupmod/proc/pupmod_src_fcd.mat','h')

s_res = squeeze(nanmean(s_res(:,:,SUBJLIST,:,:,:),6));
s_cnt = squeeze(nanmean(s_cnt(:,:,SUBJLIST,:,:,:),6));


%%
clear tperm_cnt1_n tperm_cnt1_p tperm_cnt2_n tperm_cnt2_p
clear tperm_res1_n tperm_res1_p tperm_res2_n tperm_res2_p

nperm = 10000;
v = 10;
par.subs = 200;
par.allperms = nperm/par.subs;
alpha = 0.05;
alp = 0.2;

for iperm = 1 : par.allperms
  
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,10000,v),'par')

  tperm_cnt1_n((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_cnt1_n;
  tperm_cnt1_p((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_cnt1_p;
    
  tperm_res1_n((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_res1_n;
  tperm_res1_p((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_res1_p;
    
  tperm_cnt2_n((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_cnt2_n;
  tperm_cnt2_p((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_cnt2_p;

  tperm_res2_n((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_res2_n;
  tperm_res2_p((iperm-1)*par.subs+1:(iperm)*par.subs,:)=par.tperm_res2_p;

end
%%
  for ifoi = 1 : 13

 
    [tmp1,~,~,s] = ttest(s_cnt(:,:,:,2,ifoi),s_cnt(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp1 = tmp1.*sign(s.tstat);
    [tmp2,~,~,s] = ttest(s_cnt(:,:,:,3,ifoi),s_cnt(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp2 = tmp2.*sign(s.tstat);
    [tmp3,~,~,s] = ttest(s_res(:,:,:,2,ifoi),s_res(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp3 = tmp3.*sign(s.tstat);
    [tmp4,~,~,s] = ttest(s_res(:,:,:,3,ifoi),s_res(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp4 = tmp4.*sign(s.tstat);

    temp_cnt1_n(ifoi)=sum(sum(triu(tmp1<0,1)));
    temp_cnt1_p(ifoi)=sum(sum(triu(tmp1>0,1)));
    
    temp_cnt2_n(ifoi)=sum(sum(triu(tmp2<0,1)));
    temp_cnt2_p(ifoi)=sum(sum(triu(tmp2>0,1)));

    temp_res1_n(ifoi)=sum(sum(triu(tmp3<0,1)));
    temp_res1_p(ifoi)=sum(sum(triu(tmp3>0,1)));
    
    temp_res2_n(ifoi)=sum(sum(triu(tmp4<0,1)));
    temp_res2_p(ifoi)=sum(sum(triu(tmp4>0,1)));
    
  end
  %%

% ULTIPLE COMPARISONS CORRECTION (HAWELLEK ET AL., xxxx)

for ifreq = 1 : 13
  [idx_D_cnt1_n(:,ifreq)] = floor(tiedrank(tperm_cnt1_n(:,ifreq)));
  [idx_D_cnt1_p(:,ifreq)] = floor(tiedrank(tperm_cnt1_p(:,ifreq)));
  [idx_D_cnt2_n(:,ifreq)] = floor(tiedrank(tperm_cnt2_n(:,ifreq)));
  [idx_D_cnt2_p(:,ifreq)] = floor(tiedrank(tperm_cnt2_p(:,ifreq)));
  [idx_D_res1_n(:,ifreq)] = floor(tiedrank(tperm_res1_n(:,ifreq)));
  [idx_D_res1_p(:,ifreq)] = floor(tiedrank(tperm_res1_p(:,ifreq)));
  [idx_D_res2_n(:,ifreq)] = floor(tiedrank(tperm_res2_n(:,ifreq)));
  [idx_D_res2_p(:,ifreq)] = floor(tiedrank(tperm_res2_p(:,ifreq)));
end

idx_R     = min(tperm_res2_n,[],2);
idx_D_max = sort(tperm_res2_n,'descend');
idx_D_max = idx_D_max(idx_R,:);
  

  % ----------------------------------
  % PLOT CORRECTED P-VALUES
  % ----------------------------------
  % even freqs: HC > MS
  for ifoi = 1:13
    
    corp_res2_n(ifoi)  = sum(temp_res2_n(ifoi)>idx_D_max(:,ifoi))/10000;
    corp_cnt2_n(ifoi)  = sum(temp_cnt2_n(ifoi)>idx_D_max(:,ifoi))/10000;
    
  end

%%


  %
for ifoi = 1 : 13
  
  % COMPUTE EFFECT OF DRUGS ON TASK
  p_cnt1_n(ifoi) = 1-[sum(temp_cnt1_n(ifoi)>tperm_cnt1_n(:,ifoi))]/nperm;
  p_cnt1_p(ifoi) = 1-[sum(temp_cnt1_p(ifoi)>tperm_cnt1_p(:,ifoi))]/nperm;
% 
  p_cnt2_n(ifoi) = 1-[sum(temp_cnt2_n(ifoi)>tperm_cnt2_n(:,ifoi))]/nperm;
  p_cnt2_p(ifoi) = 1-[sum(temp_cnt2_p(ifoi)>tperm_cnt2_p(:,ifoi))]/nperm;
  
  % COMPUTE EFFECT OF DRUGS ON REST
  p_res1_n(ifoi) = 1-[sum(temp_res1_n(ifoi)>tperm_res1_n(:,ifoi))]/nperm;
  p_res1_p(ifoi) = 1-[sum(temp_res1_p(ifoi)>tperm_res1_p(:,ifoi))]/nperm;

  p_res2_n(ifoi) = 1-[sum(temp_res2_n(ifoi)>tperm_res2_n(:,ifoi))]/nperm;
  p_res2_p(ifoi) = 1-[sum(temp_res2_p(ifoi)>tperm_res2_p(:,ifoi))]/nperm;
  
  
  % COMPUTE CHANGE OF CHANGE
  p_d1_p(ifoi) = 1-sum(abs(temp_res1_p(ifoi)-temp_cnt1_p(ifoi))>abs((tperm_res1_p(:,ifoi)-tperm_cnt1_p(:,ifoi))))/nperm;
  p_d1_n(ifoi) = 1-sum(abs(temp_res1_n(ifoi)-temp_cnt1_n(ifoi))>abs((tperm_res1_n(:,ifoi)-tperm_cnt1_n(:,ifoi))))/nperm;

  p_d2_p(ifoi) = 1-sum(abs(temp_res2_p(ifoi)-temp_cnt2_p(ifoi))>abs((tperm_res2_p(:,ifoi)-tperm_cnt2_p(:,ifoi))))/nperm;
  p_d2_n(ifoi) = 1-sum(abs(temp_res2_n(ifoi)-temp_cnt2_n(ifoi))>abs((tperm_res2_n(:,ifoi)-tperm_cnt2_n(:,ifoi))))/nperm;

  d1_p(ifoi) = temp_res1_p(ifoi)-temp_cnt1_p(ifoi);
  d1_n(ifoi) = temp_res1_n(ifoi)-temp_cnt1_n(ifoi);
  
  d2_p(ifoi) = temp_res2_p(ifoi)-temp_cnt2_p(ifoi);  
  d2_n(ifoi) = temp_res2_n(ifoi)-temp_cnt2_n(ifoi);
  
end

p_cnt1_n(p_cnt1_n==0) = 1-eps;
p_cnt1_p(p_cnt1_p==0) = 1-eps;
p_res1_n(p_res1_n==0) = 1-eps;
p_res1_p(p_res1_p==0) = 1-eps;

p_cnt2_n(p_cnt2_n==0) = 1-eps;
p_cnt2_p(p_cnt2_p==0) = 1-eps;
p_res2_n(p_res2_n==0) = 1-eps;
p_res2_p(p_res2_p==0) = 1-eps;

%%
p_d1_p(p_d1_p>=alpha) = nan; p_d1_p(p_d1_p<alpha) = 1; 
p_d1_n(p_d1_n>=alpha) = nan; p_d1_n(p_d1_n<alpha) = 1; 

p_d2_p(p_d2_p>=alpha) = nan; p_d2_p(p_d2_p<alpha) = 1; 
p_d2_n(p_d2_n>=alpha) = nan; p_d2_n(p_d2_n<alpha) = 1; 

%%
figure_white(12,15);

subplot(3,2,1); hold on; title('REST (Atx-Pbo)')
plot(temp_res1_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(find(p_res1_n<0.05),temp_res1_p(p_res1_n<0.05),'.','color',[0 0 0],'markersize',30)
plot(temp_res1_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(find(p_res1_p<0.05),temp_res1_p(p_res1_p<0.05),'.','color',[0 0 0],'markersize',30)

plot(1:13,repmat(1.3,[1 13]),'k-');
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]'); ylabel('# altered conn.');
set(gca,'tickdir','out'); axis([0 14 -200 1800]); axis square;
set(gca,'ytick',[0 500 1000 1500],'yticklabel',[0 500 1000 1500])

subplot(3,2,2); hold on; title('REST (Dpz-Pbo)')
plot(temp_res2_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(find(p_res2_n<0.05),temp_res2_n(p_res2_n<0.05),'.','color',[0 0 0],'markersize',30)
plot(temp_res2_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(find(p_res2_p<0.05),temp_res2_p(p_res2_p<0.05),'.','color',[0 0 0],'markersize',30)

plot(1:13,repmat(1.3,[1 13]),'k-');
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]'); ylabel('# altered conn.');
set(gca,'tickdir','out'); axis([0 14 -200 1800]); axis square;
set(gca,'ytick',[0 500 1000 1500],'yticklabel',[0 500 1000 1500])

subplot(3,2,3); hold on; title('TASK (Atx-Pbo)')

plot(temp_cnt1_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(find(p_cnt1_n<0.05),temp_cnt1_n(p_cnt1_n<0.05),'.','color',[0 0 0],'markersize',30)
plot(temp_cnt1_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(find(p_cnt1_p<0.05),temp_cnt1_p(p_cnt1_p<0.05),'.','color',[0 0 0],'markersize',30)

plot(1:13,repmat(1.3,[1 13]),'k-');
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]'); ylabel('# altered conn.');
set(gca,'tickdir','out'); axis([0 14 -200 1800]); axis square;
set(gca,'ytick',[0 500 1000 1500],'yticklabel',[0 500 1000 1500])


subplot(3,2,4); hold on; title('TASK (Dpz-Pbo)')
plot(temp_cnt2_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(find(p_cnt2_n<0.05),temp_cnt2_n(p_cnt2_n<0.05),'.','color',[0 0 0],'markersize',30)
plot(temp_cnt2_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(find(p_cnt2_p<0.05),temp_cnt2_p(p_cnt2_p<0.05),'.','color',[0 0 0],'markersize',30)

plot(1:13,repmat(1.3,[1 13]),'k-');
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]'); ylabel('# altered conn.');
set(gca,'tickdir','out'); axis([0 14 -200 1800]); axis square;
set(gca,'ytick',[0 500 1000 1500],'yticklabel',[0 500 1000 1500])
set(gca,'tickdir','out')

% NOW COMPARE ALSO THE COUNTS BETWEEN REST AND TASK

subplot(3,2,5); hold on; title('R vs T (Atx-Pbo)')

plot(d1_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(d1_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(d1_n.*p_d1_n,'.','color',[0 0 0],'markersize',30)
plot(d1_p.*p_d1_p,'.','color',[0 0 0],'markersize',30)

axis([0 14 -1000 1000])
plot(1:13,zeros([1 13]),'k-');
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]'); ylabel('Diff(# alt. conn.)');
set(gca,'tickdir','out'); axis square;

subplot(3,2,6); hold on; title('R vs T (Dpz-Pbo)')

plot(d2_n,'color',[236/255 0 139/255],'linewidth',2,'linestyle','-')
plot(d2_p,'color',[236/255 166/255 39/255],'linewidth',2,'linestyle','-')
plot(d2_n.*p_d2_n,'.','color',[0 0 0],'markersize',30)
plot(d2_p.*p_d2_p,'.','color',[0 0 0],'markersize',30)

plot(1:13,zeros([1 13]),'k-');
%
box on; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]');  ylabel('Diff(# alt. conn.)');
set(gca,'tickdir','out'); axis square;
axis([0 14 -500 500])


print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp.eps'));

%% PLOT AVERAGE FUNCTIONAL CONNECTIVITY 

for ifoi = 1 : 13
  
  [~,p_cnt1(ifoi)] = ttest(squeeze(nanmean(nanmean(s_cnt(:,:,:,2,ifoi)-s_cnt(:,:,:,1,ifoi),1),2)));
  m_cnt1(ifoi) = mean(nonzeros(triu(nanmean(s_cnt(:,:,:,2,ifoi)-s_cnt(:,:,:,1,ifoi),3),1)));
%   m_cnt1(ifoi) = mean(nonzeros(triu(nanmean(s_cnt(:,:,:,2,ifoi),3),1)))./mean(nonzeros(triu(nanmean(s_cnt(:,:,:,1,ifoi),3),1)))

  [~,p_cnt2(ifoi)] = ttest(squeeze(nanmean(nanmean(s_cnt(:,:,:,3,ifoi)-s_cnt(:,:,:,1,ifoi),1),2)));
    m_cnt2(ifoi) = mean(nonzeros(triu(nanmean(s_cnt(:,:,:,3,ifoi)-s_cnt(:,:,:,1,ifoi),3),1)));

%   m_cnt2(ifoi) = mean(nonzeros(triu(nanmean(s_cnt(:,:,:,3,ifoi),3),1)))./mean(nonzeros(triu(nanmean(s_cnt(:,:,:,1,ifoi),3),1)))
  
  s_cnt1(:,ifoi) = squeeze(nanmean(nanmean(s_cnt(:,:,:,2,ifoi)-s_cnt(:,:,:,1,ifoi),1),2));
  s_cnt2(:,ifoi) = squeeze(nanmean(nanmean(s_cnt(:,:,:,3,ifoi)-s_cnt(:,:,:,1,ifoi),1),2));
  
  [~,p_res1(ifoi)] = ttest(squeeze(nanmean(nanmean(s_res(:,:,:,2,ifoi)-s_res(:,:,:,1,ifoi),1),2)));
  m_res1(ifoi) = mean(nonzeros(triu(nanmean(s_res(:,:,:,2,ifoi)-s_res(:,:,:,1,ifoi),3),1)));
%   m_res1(ifoi) = mean(nonzeros(triu(nanmean(s_res(:,:,:,2,ifoi),3),1)))./mean(nonzeros(triu(nanmean(s_res(:,:,:,1,ifoi),3),1)))
 
  [~,p_res2(ifoi)] = ttest(squeeze(nanmean(nanmean(s_res(:,:,:,3,ifoi)-s_res(:,:,:,1,ifoi),1),2)));
  m_res2(ifoi) = mean(nonzeros(triu(nanmean(s_res(:,:,:,3,ifoi)-s_res(:,:,:,1,ifoi),3),1)));
%   m_res2(ifoi) = mean(nonzeros(triu(nanmean(s_res(:,:,:,3,ifoi),3),1)))./mean(nonzeros(triu(nanmean(s_res(:,:,:,1,ifoi),3),1)))
  
  s_res1(:,ifoi) = squeeze(nanmean(nanmean(s_res(:,:,:,2,ifoi)-s_res(:,:,:,1,ifoi),1),2));
  s_res2(:,ifoi) = squeeze(nanmean(nanmean(s_res(:,:,:,3,ifoi)-s_res(:,:,:,1,ifoi),1),2));
  
  [~, p1(ifoi)] = ttest(atanh(s_cnt1(:,ifoi)),atanh(s_res1(:,ifoi)));
  [~, p2(ifoi)] = ttest(atanh(s_cnt2(:,ifoi)),atanh(s_res2(:,ifoi)));

  
end
%
figure_white(12,8); hold on
plot(0:14,ones(15,1),'color','k','linewidth',0.5,'linestyle','-')
% plot(0:14,zeros(15,1),'color','k','linewidth',0.5,'linestyle','-')

plot(m_cnt1,'color',[1 0.2 0.1],'linewidth',3)
plot(find(p_cnt1<0.05),m_cnt1(p_cnt1<0.05),'.','color',[0 0 0],'markersize',30)
plot(m_cnt2,'color',[0.2 0.5 1],'linewidth',3)
plot(find(p_cnt2<0.05),m_cnt2(p_cnt2<0.05),'.','color',[0 0 0],'markersize',30)

plot(m_res1,'color',[1 0.2 0.1],'linewidth',3,'linestyle',':')
plot(find(p_res1<0.05),m_res1(p_res1<0.05),'.','color',[0 0 0],'markersize',30)
plot(m_res2,'color',[0.2 0.5 1],'linewidth',3,'linestyle',':')
plot(find(p_res2<0.05),m_res2(p_res2<0.05),'.','color',[0 0 0],'markersize',30)

% axis([-0 14 .65 1.35])
axis([-0 14 -0.01 0.02])

box off; set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',foi_range([1:2:13]));  xlabel('Carrier frequency [Hz]');  ylabel('Diff (Drug vs. Placebo)');
set(gca,'tickdir','out');
% legend(['ATX (Task)';'DPZ (Task)';'ATX (Rest)';'DPZ (Rest)'],'location','southeast')

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp_avg.eps'));
%% PLOT
% s = st;
% 
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);
% cmap = jet;
clim = [-0.025 0.025]
nfoi=3;
cnt = 0;
abscorr = 1;

for ifoi = 6:6
  
cnt = cnt + 1;

ss_res = s_res(:,:,:,:,ifoi);
ss_cnt = s_cnt(:,:,:,:,ifoi);

% sss_res(ifoi) = nanmean(nonzeros(triu(ss_res,1)));
% sss_(ifoi) = nanmean(nonzeros(triu(ss_res,1)));

if abscorr

  d1_res=nanmean(abs(ss_res(:,:,:,2)),3)-nanmean(abs(ss_res(:,:,:,1)),3);
  d2_res=nanmean(abs(ss_res(:,:,:,3)),3)-nanmean(abs(ss_res(:,:,:,1)),3);
  t1_res = ttest(abs(ss_res(:,:,:,2)),abs(ss_res(:,:,:,1)),'dim',3);
  t2_res = ttest(abs(ss_res(:,:,:,3)),abs(ss_res(:,:,:,1)),'dim',3);

  d1_cnt=nanmean(abs(ss_cnt(:,:,:,2)),3)-nanmean(abs(ss_cnt(:,:,:,1)),3);
  d2_cnt=nanmean(abs(ss_cnt(:,:,:,3)),3)-nanmean(abs(ss_cnt(:,:,:,1)),3);
  t1_cnt = ttest(abs(ss_cnt(:,:,:,2)),abs(ss_cnt(:,:,:,1)),'dim',3);
  t2_cnt = ttest(abs(ss_cnt(:,:,:,3)),abs(ss_cnt(:,:,:,1)),'dim',3);

else
  
  d1_res=nanmean(ss_res(:,:,:,2),3)-nanmean(ss_res(:,:,:,1),3);
  d2_res=nanmean(ss_res(:,:,:,3),3)-nanmean(ss_res(:,:,:,1),3);

  d1_cnt=nanmean(ss_cnt(:,:,:,2),3)-nanmean(ss_cnt(:,:,:,1),3);
  d2_cnt=nanmean(ss_cnt(:,:,:,3),3)-nanmean(ss_cnt(:,:,:,1),3);

end

% [t1_res p1_res]=ttest(ss_res(:,:,:,2),ss_res(:,:,:,1),'dim',3); t1_res(isnan(t1_res))=0;
% [t2_res p2_res]=ttest(ss_res(:,:,:,3),ss_res(:,:,:,1),'dim',3); t2_res(isnan(t1_res))=0;

% [t1_cnt p1_cnt]=ttest(ss_res(:,:,:,2),ss_res(:,:,:,1),'dim',3); t1_cnt(isnan(t1_cnt))=0;
% [t2_cnt p2_cnt]=ttest(ss_res(:,:,:,3),ss_res(:,:,:,1),'dim',3); t2_cnt(isnan(t1_cnt))=0;

ss_clim = squeeze(nanmean(ss_res,3)); ss_clim(ss_clim==inf)=nan;


% clim = [-0.02 0.02]
% cmap = jet;
% cmap(31:33,:)=[.95 .95 .95; .95 .95 .95;.95 .95 .95]
figure_white(15,15);

% subplot(4,3,1)
% imagesc(nanmean(ss_res(:,:,:,1),3),clim); axis square; colormap(cmap)

% subplot(4,3,2)
% imagesc(nanmean(ss_res(:,:,:,2),3),clim); axis square; 

% subplot(4,3,3)
% imagesc(nanmean(ss_res(:,:,:,3),3),clim); axis square;

subplot(1,2,1)
% imagesc(-log10(p1_res),[-2 2]); axis square; colormap(hot)
atx = tril(t1_res.*d1_res,-1);
dpz = rot90(rot90(tril(t2_res.*d2_res,-1)));
imagesc(atx+dpz,[clim]); axis square; colormap(cmap); axis off
% subplot(2,2,2)
% % imagesc(-log10(p2_res),[-2 2]); axis square; colormap(hot)
% imagesc(tril(d2_res,-1),[clim]); axis square;  axis off

% subplot(4,2,3)
% [t,p]=ttest(ss_res(:,:,:,1),ss_res(:,:,:,2),'dim',3);
% imagesc(fdr(p,0.05),[-1 1]); axis square
% 
% subplot(4,2,4)
% [t,p]=ttest(ss_res(:,:,:,1),ss_res(:,:,:,3),'dim',3);
% imagesc(fdr(tril(p),0.05),[-1 1]); axis square
% % 
% subplot(4,3,7)
% imagesc(nanmean(ss_cnt(:,:,:,1),3),clim); axis square; colormap(cmap)
% 
% subplot(4,3,8)
% imagesc(nanmean(ss_cnt(:,:,:,2),3),clim); axis square; 
% 
% subplot(4,3,9)
% imagesc(nanmean(ss_cnt(:,:,:,3),3),clim); axis square;

subplot(1,2,2)
atx = tril(t1_cnt.*d1_cnt,-1);
dpz = rot90(rot90(tril(t2_cnt.*d2_cnt,-1)));
% imagesc(-log10(p1_cnt),[-2 2]); axis square
imagesc(atx+dpz,[clim]); axis square; colormap(cmap); axis off

% imagesc(tril(d1_cnt,-1),[clim]); axis square;  axis off
% subplot(2,2,4)
% imagesc(-log10(p2_cnt),[-2 2]); axis square
% imagesc(tril(d2_cnt,-1),[clim]); axis square;  axis off

% subplot(4,2,7)
% [t,p]=ttest(ss_cnt(:,:,:,1),ss_cnt(:,:,:,2),'dim',3);
% imagesc(fdr(p,0.05),[-1 1]); axis square

% subplot(4,2,8)
% [t,p]=ttest(ss_cnt(:,:,:,1),ss_cnt(:,:,:,3),'dim',3);
% imagesc(fdr(p,0.05),[-1 1]); axis square

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_powcorr_raw_f%d_v%d.eps',ifoi,v));



end


%% PLOT THE RAW CONNECTIVITY MATRICES ACROSS CONDITIONS
figure_white(15,15)
for ifoi = 1 : 13
 subplot(4,4,ifoi)
 imagesc(nanmean(nanmean(s_res(:,:,:,ifoi),3),4),[0 0.1]); colormap('jet')
end
figure_white(15,15)
for ifoi = 1 : 13
 subplot(4,4,ifoi)
 imagesc(nanmean(nanmean(s_cnt(:,:,:,ifoi),3),4),[0 0.1]); colormap('jet')
end


%% COMPARE REST VS TASK


figure_white(20,20); hold on

  for ifoi = 1 : 13
    subplot(4,4,ifoi)
    
%     t = ttest(s_cnt(:,:,:,1,ifoi),s_res(:,:,:,1,ifoi),'dim',3);
    
    imagesc(nanmean(s_cnt(:,:,:,3,ifoi)-s_res(:,:,:,3,ifoi),3),[-0.02 0.02]); colormap(cmap)
    title(sprintf('Freq: %d Hz',foi_range(ifoi)))
    [tmp1,~,~,s] = ttest(s_cnt(:,:,:,1,ifoi),s_res(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp1 = tmp1.*sign(s.tstat);
   

    temp_n(ifoi)=sum(sum(triu(tmp1<0,1)));
    temp_p(ifoi)=sum(sum(triu(tmp1>0,1)));
    
  
    
  end
 %% 
 viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
 
 for ifoi = 1 : 13
   ifoi
   
    z   = zeros(2113,1);
    t   = ttest(s_cnt(:,:,:,1,ifoi),s_res(:,:,:,1,ifoi),'dim',3);
    par = nanmean(s_cnt(:,:,:,1,ifoi)-s_res(:,:,:,1,ifoi),3).*t;
    
   for i = 1 : 90
%      z(aalgrid.mask==i)=nanmean(nonzeros(par(i,:)));
          z(aalgrid.mask==i)=nanmean(nanmean(s_cnt(i,:,:,1,ifoi)-s_res(i,:,:,1,ifoi),3),2);

   end
   
   z(isnan(z))=0;

%    z(z~=0) = sign(z(z~=0)).*log10(abs(z(z~=0)));
   
   load sa_meg_template

   dd         = 0.75;
   g1         = sa_meg_template.grid_coarse;
   g2         = sa_meg_template.cortex10K.vc;
   par_interp = spatfiltergauss(z,g1,dd,g2);

   para = [];
   para.colorlimits = [-0.015 0.015];
    
   figure_white(10,15)
   
   for iplot = 1 : 6

     subplot(3,2,iplot)

     para.myviewdir = viewdir(iplot,:);
     a = sa_meg_template.cortex10K;

     if iplot == 5
       a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
     elseif iplot == 6
       a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
     end

     pconn_showsurface(a,para,par_interp)
     colormap(cmap)
     camlight headlight

   end
   print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp_surface_f%d_v%d.jpeg',ifoi,v));
 end

  %%
figure_white(20,20); hold on

  for ifoi = 1 : 13
    subplot(4,4,ifoi)
    
%     t = ttest(s_cnt(:,:,:,1,ifoi),s_res(:,:,:,1,ifoi),'dim',3);
    
    imagesc(nanmean(s_cnt(:,:,:,3,ifoi)-s_res(:,:,:,3,ifoi),3),[-0.02 0.02]); colormap(cmap)
    title(sprintf('Freq: %d Hz',foi_range(ifoi)))
    [tmp1,~,~,s] = ttest(s_cnt(:,:,:,1,ifoi),s_res(:,:,:,1,ifoi),'dim',3,'alpha',alp);
    tmp1 = tmp1.*sign(s.tstat);
   

    temp_n(ifoi)=sum(sum(triu(tmp1<0,1)));
    temp_p(ifoi)=sum(sum(triu(tmp1>0,1)));
    
  
    
  end
  
 print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrest_comp_pbo.pdf'));
  
 figure_white(20,20); hold on
 plot(temp_n)
 plot(temp_p)

%%
figure_white(25,20); 
for ifoi  = 1 : 13
      subplot(3,5,ifoi); hold on

for m = 1:2
    if m == 1
      s1 = nonzeros(triu(nanmean(abs(s_res(:,:,:,1,ifoi),3),1));
      s2 = nonzeros(triu(nanmean(s_res(:,:,:,2,ifoi),3),1));
    else
      s1 = nonzeros(triu(nanmean(s_cnt(:,:,:,1,ifoi),3),1));
      s2 = nonzeros(triu(nanmean(s_cnt(:,:,:,2,ifoi),3),1));
    end
    corr(s1,s2)
    scatter(s1,s2,'.'); axis square; box on; 
    set(gca,'tickdir','out');
    axis([-0.01 0.14 -0.01 0.14])

  end
end