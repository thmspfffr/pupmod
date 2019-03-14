%% PLOT NUMBER OF ALTERED CORRELATIONS, INCLUDING STATISTICS
% --------------------------
% This script obtains the empircal number of altered correlations, by
% calling pupmod_compute_altered_correlations.m and obtains a corrected 
% p-values from a permutation distribution (computed in
% pupmod_src_powcorr_permtest.m). The actual p-values are obtained calling
% the function pupmod_all_powcorr_getstatistics.m)
% --------------------------
% CONTENTS
% --------------
% (1) PLOT: No stats
% (2) PLOT: P-Values (corrected) 
% (3) PLOT: Altered correlations
% (4) Circular graphs
% (5) Plot altered correlations (per voxel)
% --------------

clear

% -------------
% version of cleaned data: 
% v1: AAL, v2: 400 vertices (cortex)
% -------------
v = 12;
% -------------

% load  data
cd ~/pupmod/matlab/
para.avg = 0;
% load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
cleandat = pupmod_loadpowcorr(v,para);
%%

para = [];
para.nfreq = 1:13;
para.alpha = 0.05;

emp = pupmod_compute_altered_correlations(cleandat,para);

%% BAR PLOT AT F = 6 and F = 7
figure; set(gcf,'color','w');
subplot(2,2,3); hold on
bar([1 2],[emp.n_p_atx(6,1)-emp.n_p_atx(6,2) emp.n_n_atx(6,1)-emp.n_n_atx(6,2)]); axis square
axis(gca,[0.5 2.5 -0.6 0.6]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots

subplot(2,2,4); hold on
bar([1 2],[emp.n_p_dpz(7,1)-emp.n_p_dpz(7,2) emp.n_n_dpz(7,1)-emp.n_n_dpz(7,2)]); axis square
axis(gca,[0.5 2.5 -0.6 0.6]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_alteredcorr_bar_v%d.eps',v))

% CHANGE in CORRELATION (PERCENT)
fc= squeeze(nanmean(nanmean(nanmean(cleandat),2),3));

atx_res = 100*(fc(2,1,6)-fc(1,1,6))/fc(1,1,6);
atx_cnt = 100*(fc(2,2,6)-fc(1,2,6))/fc(1,2,6);
dpz_res = 100*(fc(3,1,7)-fc(1,1,7))/fc(1,1,7);
dpz_cnt = 100*(fc(3,2,7)-fc(1,2,7))/fc(1,2,7);
atx_ctx = atx_res - atx_cnt;
dpz_ctx = dpz_res - dpz_cnt;

fc = squeeze(nanmean(nanmean(cleandat,2)));

s_atx_res = 100*std((fc(:,2,1,6)-fc(:,1,1,6))./fc(:,1,1,6),1)/sqrt(size(fc,1))
s_atx_cnt = 100*std((fc(:,2,2,6)-fc(:,1,2,6))./fc(:,1,2,6),1)/sqrt(size(fc,1))
s_dpz_res = 100*std((fc(:,3,1,7)-fc(:,1,1,7))./fc(:,1,1,7),1)/sqrt(size(fc,1))
s_dpz_cnt = 100*std((fc(:,3,2,7)-fc(:,1,2,7))./fc(:,1,2,7),1)/sqrt(size(fc,1))
s_atx_ctx = std(100*(((fc(:,2,1,6)-fc(:,1,1,6))./fc(:,1,1,6))-((fc(:,2,2,6)-fc(:,1,2,6))./fc(:,1,2,6))))/sqrt(size(fc,1))
s_dpz_ctx = std(100*(((fc(:,3,1,7)-fc(:,1,1,7))./fc(:,1,1,7))-((fc(:,3,2,7)-fc(:,1,2,7))./fc(:,1,2,7))))/sqrt(size(fc,1))

figure; set(gcf,'color','w');
subplot(2,2,3); hold on
bar([1 2 3],[atx_res atx_cnt atx_ctx]); axis square
line([1 1],[atx_res-s_atx_res atx_res+s_atx_res],'color','k')
line([2 2],[atx_cnt-s_atx_cnt atx_cnt+s_atx_cnt],'color','k')
line([3 3],[atx_ctx-s_atx_ctx atx_ctx+s_atx_ctx],'color','k')
axis(gca,[0.5 3.5 -30 30]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots

subplot(2,2,4); hold on
bar([1 2 3],[dpz_res dpz_cnt dpz_ctx]); axis square
line([1 1],[dpz_res-s_dpz_res dpz_res+s_dpz_res],'color','k')
line([2 2],[dpz_cnt-s_dpz_cnt dpz_cnt+s_dpz_cnt],'color','k')
line([3 3],[dpz_ctx-s_dpz_ctx dpz_ctx+s_dpz_ctx],'color','k')
axis(gca,[0.5 3.5 -30 30]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_rawcorr_bar_v%d.eps',v))



%% (1) PLOT: No stats
% ------------------

linewidth = 2;

figure; set(gcf,'color','w');

subplot(3,2,1); hold on

plot(emp.n_p_atx(:,1),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
ylabel('Frac. of altered corr. [%]')

subplot(3,2,3); hold on

plot(emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
ylabel('Frac. of altered corr. [%]')

subplot(3,2,2); hold on

plot(emp.n_p_dpz(:,1),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,4); hold on

plot(emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,5); hold on

plot(emp.n_p_atx(:,1)-emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1)-emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')

subplot(3,2,6); hold on

plot(emp.n_p_dpz(:,1)-emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1)-emp.n_p_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')

%%
if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
  % Settings:
  para = [];
  para.nfreq = [1:13]; % freqs 1 - 13
  para.alpha = 0.05; % alpha for adjacency
  para.ver = v;
  para.nperm = 10000;
  para.nsubs = 250;
  para.type = 'global';
  para.cond = 'atx';
  para.allperms = para.nperm/para.nsubs;

  [outp_atx, emp_atx] = pupmod_src_powcorr_getstatistics(para);
  para.cond = 'dpz';
  [outp_dpz, emp_dpz] = pupmod_src_powcorr_getstatistics(para);
  
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v),'outp_atx','outp_dpz','emp_atx','emp_dpz')

%   for iperm = 1 : para.allperms
%   
%     fprintf('Load permutation distributions: %d / %d ...\n',iperm,para.allperms)
%   
%     load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,para.nperm,para.ver),'par')
%     
%     atx(:,(iperm-1)*para.nsubs+1:(iperm)*para.nsubs) = par.tperm_cnt1_pervoxel_p(:,:,6);
%     p_atx_vox = 1-sum(abs(emp.n_p_atx_pervoxel(:,6,2))>abs(atx),2)/para.nperm;
%     
%     dpz(:,(iperm-1)*para.nsubs+1:(iperm)*para.nsubs) = par.tperm_res2_pervoxel_p(:,:,7);
%     p_dpz_vox = 1-sum(abs(emp.n_n_dpz_pervoxel(:,7,1))>abs(dpz),2)/para.nperm;
%     
%   ends
  



else
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
end
%% (2) PLOT: P-Values (corrected) 

figure;  set(gcf,'color','white')

subplot(3,2,1); hold on
plot(-log10(outp_atx.p_res1_n),'b-','linewidth',3)
plot(-log10(outp_atx.p_res1_p),'r-','linewidth',3)

line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Rest')

subplot(3,2,2); hold on
plot(-log10(outp_dpz.p_res2_n),'b-','linewidth',3)
plot(-log10(outp_dpz.p_res2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

subplot(3,2,3); hold on
plot(-log10(outp_atx.p_cnt1_n),'b-','linewidth',3)
plot(-log10(outp_atx.p_cnt1_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Task')

subplot(3,2,4); hold on
plot(-log10(outp_dpz.p_cnt2_n),'b-','linewidth',3)
plot(-log10(outp_dpz.p_cnt2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

% subplot(3,2,5); hold on
% plot(-log10(p_d1_n),'b-','linewidth',3)
% plot(-log10(p_d1_p),'r-','linewidth',3)
% line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
% axis([0 14 0 4])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
% 
% subplot(3,2,6); hold on
% plot(-log10(p_d2_n),'b-','linewidth',3)
% plot(-log10(p_d2_p),'r-','linewidth',3)
% line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
% axis([0 14 0 4])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
% % tp_editplots
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp_pval.pdf'));

%% (3a) PLOT: Altered correlations
% Plot altered correlations and highlights significant differences
% --------------
if v == 1
  lims = [0 0.25 0.5 0.75];
  lims_lab = num2cell([0 25 50 75]);
  lims_ctx = [-1 -0.5 0; 0 0.5 1];
else
  lims = [0 0.25 0.5];
  lims_lab = num2cell([0 25 50]);
  lims_ctx = [-0.5 0 0.5; -0.5 0 0.5];
end
  
alpha1 = 0.01;
alpha2 = 0.01;
alpha3 = 0.001;

figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(emp_atx.n_p_atx(:,1),'r-','linewidth',2)
plot(emp_atx.n_n_atx(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 14 lims(1)-0.05 lims(end)+0.175])
plot(find(outp_atx.p_res1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
tp_editplots
pos(1,:)=get(gca,'Position')

subplot(4,2,2); hold on
plot(emp_atx.n_p_dpz(:,1),'r-','linewidth',2)
plot(emp_atx.n_n_dpz(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
axis([0 14 lims(1)-0.05 lims(end)+0.175])
plot(find(outp_atx.p_res2_p<alpha1),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_n<alpha1),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_p<alpha2),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_n<alpha2),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_p<alpha3),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res2_n<alpha3),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
tp_editplots
pos(2,:)=get(gca,'Position')

subplot(4,2,3); hold on
plot(emp_atx.n_p_atx(:,2),'r-','linewidth',2)
plot(emp_atx.n_n_atx(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 14 lims(1)-0.05 lims(end)+0.175])
plot(find(outp_atx.p_cnt1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots
pos(3,:)=get(gca,'Position')

subplot(4,2,4); hold on
plot(emp_atx.n_p_dpz(:,2),'r-','linewidth',2)
plot(emp_atx.n_n_dpz(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
axis([0 14 lims(1)-0.05 lims(end)+0.175])
plot(find(outp_atx.p_cnt2_p<alpha1),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_n<alpha1),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_p<alpha2),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_n<alpha2),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_p<alpha3),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt2_n<alpha3),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots
pos(4,:)=get(gca,'Position')

subplot(4,2,5); hold on
plot(emp_atx.n_n_context_atx,'b-','linewidth',2)
plot(emp_atx.n_p_context_atx,'r-','linewidth',2)
axis([0 14 lims_ctx(1,1)-0.1 lims_ctx(1,end)+0.1])
set(gca,'tickdir','out','ytick',[-0.5 0 0.5],'yticklabel',num2cell([-50 0 50]))
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(outp_atx.p_context1_p<alpha1),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_n<alpha1),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_p<alpha2),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_n<alpha2),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_p<alpha3),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context1_n<alpha3),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots
pos(5,:)=get(gca,'Position')

subplot(4,2,6); hold on
plot(emp_atx.n_n_context_dpz,'b-','linewidth',2)
plot(emp_atx.n_p_context_dpz,'r-','linewidth',2)
axis([0 14 lims_ctx(2,1)-0.1 lims_ctx(2,end)+0.1])
set(gca,'tickdir','out','ytick',[-0.5 0 0.5],'yticklabel',num2cell([-50 0 50]))
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
xlabel('Carrier frequency [Hz]'); 
plot(find(outp_atx.p_context2_p<alpha1),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_n<alpha1),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_p<alpha2),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_n<alpha2),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_p<alpha3),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context2_n<alpha3),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots
pos(6,:)=get(gca,'Position')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_powcorr_drugeffects_lineplots_v%d.pdf',v));

%% (5) PLOT ALTERED CORRELATIONS (PER VOXEL)

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

para          = [];
para.nfreq    = [1:13]; % freqs 1 - 13
para.alpha    = 0.05; % alpha for adjacency
para.ver      = v;
para.nperm    = 10000;
para.nsubs    = 250;
para.type     = 'local';
para.allperms = para.nperm/para.nsubs;
para.emp      = emp;
para.cleaned  = 0;
% ---------
% Obtain stats for atomoxetine condition
para.cond     = 'atx';
[outp_atx]    = pupmod_src_powcorr_getstatistics(para);
% ---------
% Obtain stats for donepezil condition
para.cond     = 'dpz';
[outp_dpz]    = pupmod_src_powcorr_getstatistics(para);
% ---------
% Obtain stats for task vs rest
para.cond     = 'taskvsrest';
[outp_tvr]    = pupmod_src_powcorr_getstatistics(para);
% ---------
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
%% PLOT RESULTS FOR ATX (TASK) AND DPZ (REST)

close all

ifoi = 6; icond = 2;

cmap = autumn;

par = emp.n_p_atx_pervoxel(:,ifoi,icond);
par(outp_atx.pval_p_atx(:,icond,ifoi)>=0.025) = 0;

cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-0.75 0.75];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_atx_f%d_c%d_v%d.png',ifoi,icond,v);
tp_plot_surface(par,sa_template,para)


ifoi = 7; icond = 1;

cmap = autumn;
cmap(:,1) = 0; cmap(:,3) = 1;

par = emp.n_n_dpz_pervoxel(:,ifoi,icond);
par(outp_dpz.pval_n_dpz(:,icond,ifoi)>=0.025) = 0;

cmap = [cmap(:,:); 0.98*ones(1,3); cmap(:,:)];
para = [];
para.clim = [-0.75 0.75];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_dpz_f%d_c%d_v%d.png',ifoi,icond,v);
tp_plot_surface(par,sa_template,para)

% TASK VS REST
% ----------------------

ifoi = 6; icond = 1;

cmap = autumn;
cmap(:,1) = 0; cmap(:,3) = 1;

par = emp.taskvsrest_n_pervoxel(:,ifoi);
par(outp_tvr.pval_n_tvr(:,ifoi)>=0.025) = 0;

cmap = [cmap(:,:); 0.98*ones(1,3); cmap(:,:)];
para      = [];
para.clim = [-1 1];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_tvr_f%d_c%d_v%d.png',ifoi,icond,v);
tp_plot_surface(par,sa_template,para)



%%
fprintf('------------\n')

icond = 2;
ifoi= 6;

g2    = sa_template.cortex10K.vc(sa_template.cortex10K.vc(:,1)>0,:);
dd    = 0.5;

grid_lh = grid(grid(:,1)<0,:); grid_lh(:,1)=abs(grid_lh(:,1));
idx_lh = grid(:,1)<0;

grid_rh = grid(grid(:,1)>0,:); 
idx_rh = grid(:,1)>0;

par = emp.n_n_dpz_pervoxel(idx_rh,ifoi,icond);
par1 = spatfiltergauss(par,grid_rh,dd,g2);

par = emp.n_n_dpz_pervoxel(idx_lh,ifoi,icond);
par2 = spatfiltergauss(par,grid_lh,dd,g2);

[r,p]=corr(par1,par2);
fprintf('Corr (DPZ): r = %.3f\n', r)

par = emp.n_p_atx_pervoxel(idx_rh,ifoi,icond);
par1 = spatfiltergauss(par,grid_rh,dd,g2);

par = emp.n_p_atx_pervoxel(idx_lh,ifoi,icond);
par2 = spatfiltergauss(par,grid_lh,dd,g2);

[r,p]=corr(par1,par2);
fprintf('Corr (ATX): r = %.3f\n', r)

%% PLOT TASK VS REST

%% (3a) PLOT: Altered correlations
% Plot altered correlations and highlights significant differences
% --------------
if v == 1
  lims = [0 0.25 0.5 0.75];
  lims_lab = num2cell([0 25 50 75]);
  lims_ctx = [-1 -0.5 0; 0 0.5 1];
else
  lims = [0 0.25 0.5];
  lims_lab = num2cell([0 25 50]);
  lims_ctx = [-0.5 0 0.5; -0.5 0 0.5];
end
  
alpha1 = 0.01;
alpha2 = 0.01;
alpha3 = 0.001;

figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(emp_atx.taskvsrest_p,'r-','linewidth',2)
plot(emp_atx.taskvsrest_n,'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 14 lims(1)-0.05 lims(end)+0.175])
plot(find(outp_atx.p_tvr_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_tvr_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_tvr_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_tvr_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_tvr_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_tvr_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_tvr_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_tvr_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_tvr_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_tvr_p<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_tvr_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_tvr_n<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
tp_editplots
pos(1,:)=get(gca,'Position')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_powcorr_tvr_lineplots_v%d.pdf',v));

% SURFACE PLOTS
% ======================
% TASK VS REST
% ----------------------

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

para          = [];
para.nfreq    = [1:13]; % freqs 1 - 13
para.alpha    = 0.05; % alpha for adjacency
para.ver      = v;
para.nperm    = 10000;
para.nsubs    = 250;
para.type     = 'local';
para.allperms = para.nperm/para.nsubs;
para.emp      = emp;
para.cleaned  = 0;
% ---------
% Obtain stats for task vs rest
para.cond     = 'taskvsrest';
[outp_tvr]    = pupmod_src_powcorr_getstatistics(para);
% ---------
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/


ifoi = 2; icond = 1;

cmap = autumn;
cmap(:,1) = 0; cmap(:,3) = 1;

par = emp.taskvsrest_p_pervoxel(:,ifoi);
par(outp_tvr.pval_p_tvr(:,ifoi)>=0.05) = 0;

cmap = [cmap(:,:); 0.98*ones(1,3); cmap(:,:)];
para      = [];
para.clim = [-1 1];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_tvr_f%d_c%d_v%d.png',ifoi,icond,v);
tp_plot_surface(par,sa_template,para)

