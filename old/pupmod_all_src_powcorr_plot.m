%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 12;
ALPHA = 0.05;
para.nfreq = 13;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

ord = pconn_randomization;

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

para.fcsize = size(cleandat,1);
para.nfreq = 13;
para.alpha = 0.05;

emp = pupmod_compute_altered_correlations(cleandat,para);


clear cleandat s_fc

% --------------------------
% load stats
% --------------------------
% see pupmod_src_powcorr_permtest.m
% see pupmod_src_powcorr_getstatistics
% --------------------------

para.nsubs = 100;
para.nperm = 50000;

load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_corrected_pvalues_subs%d_nperm%d_v%d.mat',para.nsubs,para.nperm,v))

%% PLOT NUMBER OF ALTERED CORRELATIONS (across all voxels)

figure; set(gcf,'color','white')

subplot(3,2,1); hold on
plot(emp.n_p_atx(:,1),'r-','linewidth',3)
plot(emp.n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
ylabel('Altered corr. [%]')
title('Rest')
axis([0 14 -0.05 0.55])
plot(find(outp.atx.p_res1_p<0.025),emp.n_p_atx(find(outp.atx.p_res1_p<0.025),1),'k.','markersize',30)
plot(find(outp.atx.p_res1_n<0.025),emp.n_n_atx(find(outp.atx.p_res1_n<0.025),1),'k.','markersize',30)
tp_editplots
% plot(prctile(perm_n_n_atx(:,:,1),95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(perm_n_p_atx(:,:,1),95),'linewidth',1,'color','r','linestyle',':')

subplot(3,2,2); hold on
plot(emp.n_p_dpz(:,1),'r-','linewidth',3)
plot(emp.n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([0 14 -0.05 0.55])
plot(find(outp.dpz.p_res2_p<0.025),emp.n_p_dpz(find(outp.dpz.p_res2_p<0.025),1),'k.','markersize',30)
plot(find(outp.dpz.p_res2_n<0.025),emp.n_n_dpz(find(outp.dpz.p_res2_n<0.025),1),'k.','markersize',30)
tp_editplots
% plot(prctile(perm_n_n_dpz(:,:,1),95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(perm_n_p_dpz(:,:,1),95),'linewidth',1,'color','r','linestyle',':')

subplot(3,2,3); hold on
plot(emp.n_p_atx(:,2),'r-','linewidth',3)
plot(emp.n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
ylabel('Altered corr. [%]')
title('Task')
axis([0 14 -0.05 0.55])
plot(find(outp.atx.p_cnt1_p<0.025),emp.n_p_atx(find(outp.atx.p_cnt1_p<0.025),2),'k.','markersize',30)
plot(find(outp.atx.p_cnt1_n<0.025),emp.n_n_atx(find(outp.atx.p_cnt1_n<0.025),2),'k.','markersize',30)
tp_editplots
% plot(prctile(perm_n_n_atx(:,:,2),95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(perm_n_p_atx(:,:,2),95),'linewidth',1,'color','r','linestyle',':')

subplot(3,2,4); hold on
plot(emp.n_p_dpz(:,2),'r-','linewidth',3)
plot(emp.n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([0 14 -0.05 0.55])
plot(find(outp.dpz.p_cnt2_p<0.025),emp.n_p_dpz(find(outp.dpz.p_cnt2_p<0.025),2),'k.','markersize',30)
plot(find(outp.dpz.p_cnt2_n<0.025),emp.n_n_dpz(find(outp.dpz.p_cnt2_n<0.025),2),'k.','markersize',30)
tp_editplots
% plot(prctile(perm_n_n_dpz(:,:,2),95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(perm_n_p_dpz(:,:,2),95),'linewidth',1,'color','r','linestyle',':')

subplot(3,2,5); hold on
plot(emp.n_n_atx(:,2)-emp.n_n_atx(:,1),'b-','linewidth',3)
plot(emp.n_p_atx(:,2)-emp.n_p_atx(:,1),'r-','linewidth',3)
axis([0 14 -0.55 0.55])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
% plot(find(p_d1_p<0.025),d1_p(find(p_d1_p<0.025)),'k.','markersize',30)
% plot(find(p_d1_n<0.025),d1_n(find(p_d1_n<0.025)),'k.','markersize',30)
tp_editplots

subplot(3,2,6); hold on
plot(emp.n_n_dpz(:,2)-emp.n_n_dpz(:,1),'b-','linewidth',3)
plot(emp.n_p_dpz(:,2)-emp.n_p_dpz(:,1),'r-','linewidth',3)
axis([0 14 -0.55 0.55])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]');
% plot(find(p_d2_p<0.025),d2_p(find(p_d2_p<0.025)),'k.','markersize',30)
% plot(find(p_d2_n<0.025),d2_n(find(p_d2_n<0.025)),'k.','markersize',30)
tp_editplots

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_alteredcorr_global_v%d.pdf',v));

%% PLOT SPATIAL DISTRIBUTION







