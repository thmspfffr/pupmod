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
v = 3;
% -------------

% load  data
cd ~/pupmod/matlab/

% para.avg = 0;

fc = pupmod_loadpowcorr(v,SUBJLIST,1);

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_within_v%d.mat',v));
%%
para = [];
para.nfreq = 1:17;
para.alpha = 0.05;
para.avg = 0;
% fc = nanmean(fc,7);
% if para.avg == 1
emp = pupmod_compute_altered_correlations(cleandat,para);
% else
%   iblock = 1;
%   emp = pupmod_compute_altered_correlations(fc(:,:,:,:,:,:,iblock),para);
% end

%% (1) PLOT: No stats
% ------------------

linewidth = 1;

figure; set(gcf,'color','w');

subplot(4,2,1); hold on

plot(emp.n_p_atx(:,1),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.9])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))

ylabel('Frac. of altered corr. [%]')
tp_editplots

subplot(4,2,3); hold on

plot(emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.9])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))

ylabel('Frac. of altered corr. [%]')
tp_editplots

subplot(4,2,2); hold on

plot(emp.n_p_dpz(:,1),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.9])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))

tp_editplots

subplot(4,2,4); hold on

plot(emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.9])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))

tp_editplots

subplot(4,2,5); hold on

plot(emp.n_p_atx(:,1)-emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1)-emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -1 1])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-1 0 1],'yticklabel',num2cell([-1 0 1]))

xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')
tp_editplots

subplot(4,2,6); hold on

plot(emp.n_p_dpz(:,1)-emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1)-emp.n_p_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -1 1])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-1 0 1],'yticklabel',num2cell([-1 0 1]))

xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_powcorr_drugeffects_nostats_lineplots_b%d.pdf',iblock));

%%
% if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
  % Settings:
  para = [];
  para.nfreq = [1:17]; % freqs 1 - 13
  para.alpha = 0.05; % alpha for adjacency
  para.ver = v;
  para.nperm = 10000;
  para.nsubs = 200;
  para.type = 'global';
  para.cond = 'atx';
 para.correction_method = 'single_threshold';
 para.emp = emp;
  para.allperms = para.nperm/para.nsubs;

  [outp_atx] = pupmod_src_powcorr_getstatistics_cleaned(para);
%   para.cond = 'dpz';
%   [outp_dpz] = pupmod_src_powcorr_getstatistics_cleaned(para);
  
%   save(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v),'outp_atx','outp_dpz','emp_atx','emp_dpz')

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
%   end
  



% else
%   load(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
% end
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

alpha1 = 0.025;
alpha2 = 0.01;
alpha3 = 0.001;

figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(emp_atx.n_p_atx(:,1),'r-','linewidth',3)
plot(emp_atx.n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))
ylabel('Altered corr. [%]')
title('Rest')
axis([0 14 -0.05 0.9])
plot(find(outp_atx.p_res1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(4,2,2); hold on
plot(emp_atx.n_p_dpz(:,1),'r-','linewidth',3)
plot(emp_atx.n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))
axis([0 14 -0.05 0.9])
plot(find(outp_atx.p_res2_p<alpha1),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_n<alpha1),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_p<alpha2),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_n<alpha2),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_p<alpha3),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res2_n<alpha3),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(4,2,3); hold on
plot(emp_atx.n_p_atx(:,2),'r-','linewidth',3)
plot(emp_atx.n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))
ylabel('Altered corr. [%]')
title('Task')
axis([0 14 -0.05 0.9])
plot(find(outp_atx.p_cnt1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(4,2,4); hold on
plot(emp_atx.n_p_dpz(:,2),'r-','linewidth',3)
plot(emp_atx.n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.25 0.50 0.75],'yticklabel',num2cell([0 0.25 0.50 0.75]))
axis([0 14 -0.05 0.9])
plot(find(outp_atx.p_cnt2_p<alpha1),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_n<alpha1),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_p<alpha2),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_n<alpha2),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_p<alpha3),emp_atx.n_p_dpz(find(outp_dpz.p_cnt2_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt2_n<alpha3),emp_atx.n_n_dpz(find(outp_dpz.p_cnt2_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(4,2,5); hold on
plot(emp_atx.n_n_context_atx,'b-','linewidth',3)
plot(emp_atx.n_p_context_atx,'r-','linewidth',3)
axis([0 14 -1 1])
set(gca,'tickdir','out','ytick',[-1 0 1],'yticklabel',num2cell([-1 0 1]))
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(outp_atx.p_context1_p<alpha1),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_n<alpha1),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_p<alpha2),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_n<alpha2),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_p<alpha3),emp_atx.n_p_context_atx(find(outp_atx.p_context1_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context1_n<alpha3),emp_atx.n_n_context_atx(find(outp_atx.p_context1_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(4,2,6); hold on
plot(emp_atx.n_n_context_dpz,'b-','linewidth',3)
plot(emp_atx.n_p_context_dpz,'r-','linewidth',3)
axis([0 14 -1 1])
set(gca,'tickdir','out','ytick',[-1 0 1],'yticklabel',num2cell([-1 0 1]))
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
xlabel('Carrier frequency [Hz]'); 
plot(find(outp_atx.p_context2_p<alpha1),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_n<alpha1),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_p<alpha2),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_n<alpha2),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_p<alpha3),emp_atx.n_p_context_dpz(find(outp_dpz.p_context2_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context2_n<alpha3),emp_atx.n_n_context_dpz(find(outp_dpz.p_context2_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_powcorr_drugeffects_lineplots.pdf'));

%% (5) PLOT ALTERED CORRELATIONS (PER VOXEL)

figure; set(gcf,'color','w');
par = 100*(nanmean(cleandat(:,:,:,2,2,7),3)-nanmean(cleandat(:,:,:,1,2,7),3))./nanmean(cleandat(:,:,:,1,2,7),3);
par1 = triu(par,1);

par = 100*(nanmean(cleandat(:,:,:,2,1,7),3)-nanmean(cleandat(:,:,:,1,1,7),3))./nanmean(cleandat(:,:,:,1,2,7),3);

imagesc(par,[0 30])
colormap(redblue); axis off



