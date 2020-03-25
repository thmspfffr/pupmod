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
% (1) PLOT: P-Values (corrected) 
% (2) PLOT: Altered correlations
% (4) Plot altered correlations (per voxel)
% -------------

clear

% -------------
% version of cleaned data: 
% v1: cortex 400 vertices, v2: 400 vertices (cortex)
% -------------
v = 3;
% -------------
%%
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

fc = pupmod_loadpowcorr(v,SUBJLIST,1);

para = [];
para.nfreq = 1:17;
para.alpha = 0.05;

emp = pupmod_compute_altered_correlations(fc,para);

%% GET STATISTICS 

if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
  % Settings:
  para = [];
  para.nfreq = [1:size(emp.n_p_atx,1)]; % freqs 1 - 17
  para.alpha = 0.05; % alpha for adjacency
  para.ver = v;
  para.nperm = 10000;
  para.nsubs = 200;
  para.type = 'global';
  para.cond = 'atx';
  para.correction_method = 'single_threshold';
  para.emp = emp;
  
  outp_atx = pupmod_src_powcorr_getstatistics(para);
  
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v),'outp_atx','emp')

else
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
end

%% (3a) PLOT: Altered correlations
% Plot altered correlations and highlights significant differences
% --------------
nfreq = 17;
freq_start = 1;

if v == 1
  lims = [0 0.25];
  lims_lab = num2cell([0 .25]);
  lims_ctx = [-0.25 0; 0 0.25];
else
  lims = [0 0.25];
  lims_lab = num2cell([0 .25]);
  lims_ctx = [-0.25 0; 0 0.25];
end

markersize = 4;
  
alpha1 = 0.05;
alpha2 = 0.01;
alpha3 = 0.001;

figure; set(gcf,'color','w')

subplot(4,3,1); hold on
plot(freq_start:nfreq,emp.n_p_atx(freq_start:nfreq,1),'r-','linewidth',2)
plot(freq_start:nfreq,emp.n_n_atx(freq_start:nfreq,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
set(gca,'tickdir','out','ytick',cell2mat(lims_lab),'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 nfreq lims(1)-0.05 lims(end)])

plot(find(outp_atx.p_res1_p<alpha1),emp.n_p_atx(find(outp_atx.p_res1_p<alpha1),1),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_res1_n<alpha1),emp.n_n_atx(find(outp_atx.p_res1_n<alpha1),1),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_res1_p<alpha2),emp.n_p_atx(find(outp_atx.p_res1_p<alpha2),1),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_res1_n<alpha2),emp.n_n_atx(find(outp_atx.p_res1_n<alpha2),1),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_res1_p<alpha3),emp.n_p_atx(find(outp_atx.p_res1_p<alpha3),1),'ro','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_res1_n<alpha3),emp.n_n_atx(find(outp_atx.p_res1_n<alpha3),1),'ro','markersize',markersize,'markerfacecolor','m')

tp_editplots
pos(1,:)=get(gca,'Position')
axis([freq_start-1 nfreq lims(1)-0.05 lims(end)])

subplot(4,3,2); hold on
plot(freq_start:nfreq,emp.n_p_dpz(freq_start:nfreq,1),'r-','linewidth',2)
plot(freq_start:nfreq,emp.n_n_dpz(freq_start:nfreq,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
set(gca,'tickdir','out','ytick',cell2mat(lims_lab),'yticklabel',lims_lab)
axis([0 size(emp.n_p_atx,1) lims(1)-0.05 lims(end)])
plot(find(outp_atx.p_res2_p<alpha1),emp.n_p_dpz(find(outp_atx.p_res2_p<alpha1),1),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_res2_n<alpha1),emp.n_n_dpz(find(outp_atx.p_res2_n<alpha1),1),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_res2_p<alpha2),emp.n_p_dpz(find(outp_atx.p_res2_p<alpha2),1),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_res2_n<alpha2),emp.n_n_dpz(find(outp_atx.p_res2_n<alpha2),1),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_res2_p<alpha3),emp.n_p_dpz(find(outp_atx.p_res2_p<alpha3),1),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_res2_n<alpha3),emp.n_n_dpz(find(outp_atx.p_res2_n<alpha3),1),'ko','markersize',markersize,'markerfacecolor','m')
tp_editplots
pos(2,:)=get(gca,'Position')
axis([freq_start-1 nfreq lims(1)-0.05 lims(end)])

subplot(4,3,4); hold on
plot(freq_start:nfreq,emp.n_p_atx(freq_start:nfreq,2),'r-','linewidth',2)
plot(freq_start:nfreq,emp.n_n_atx(freq_start:nfreq,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
set(gca,'tickdir','out','ytick',cell2mat(lims_lab),'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 size(emp.n_p_atx,1) lims(1)-0.05 lims(end)])
plot(find(outp_atx.p_cnt1_p<alpha1),emp.n_p_atx(find(outp_atx.p_cnt1_p<alpha1),2),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_n<alpha1),emp.n_n_atx(find(outp_atx.p_cnt1_n<alpha1),2),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_p<alpha2),emp.n_p_atx(find(outp_atx.p_cnt1_p<alpha2),2),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_n<alpha2),emp.n_n_atx(find(outp_atx.p_cnt1_n<alpha2),2),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_p<alpha3),emp.n_p_atx(find(outp_atx.p_cnt1_p<alpha3),2),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_cnt1_n<alpha3),emp.n_n_atx(find(outp_atx.p_cnt1_n<alpha3),2),'ko','markersize',markersize,'markerfacecolor','m')
tp_editplots
pos(3,:)=get(gca,'Position')
axis([freq_start-1 nfreq lims(1)-0.05 lims(end)])

subplot(4,3,5); hold on
plot(freq_start:nfreq,emp.n_p_dpz(freq_start:nfreq,2),'r-','linewidth',2)
plot(freq_start:nfreq,emp.n_n_dpz(freq_start:nfreq,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
set(gca,'tickdir','out','ytick',cell2mat(lims_lab),'yticklabel',lims_lab)
axis([0 size(emp.n_p_atx,1) lims(1)-0.05 lims(end)])
plot(find(outp_atx.p_cnt2_p<alpha1),emp.n_p_dpz(find(outp_atx.p_cnt2_p<alpha1),2),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_n<alpha1),emp.n_n_dpz(find(outp_atx.p_cnt2_n<alpha1),2),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_p<alpha2),emp.n_p_dpz(find(outp_atx.p_cnt2_p<alpha2),2),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_n<alpha2),emp.n_n_dpz(find(outp_atx.p_cnt2_n<alpha2),2),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_p<alpha3),emp.n_p_dpz(find(outp_atx.p_cnt2_p<alpha3),2),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_cnt2_n<alpha3),emp.n_n_dpz(find(outp_atx.p_cnt2_n<alpha3),2),'ko','markersize',markersize,'markerfacecolor','m')
tp_editplots
pos(4,:)=get(gca,'Position')
axis([freq_start-1 nfreq lims(1)-0.05 lims(end)])

subplot(4,3,7); hold on
plot(freq_start:nfreq,emp.n_n_context_atx(freq_start:nfreq),'b-','linewidth',2)
plot(freq_start:nfreq,emp.n_p_context_atx(freq_start:nfreq),'r-','linewidth',2)
axis([freq_start-1 nfreq lims_ctx(1,1) lims_ctx(1,end)])
set(gca,'tickdir','out','ytick',[-0.25 0 0.25],'yticklabel',num2cell([-25 0 25]))
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(outp_atx.p_context_atx_p<alpha1),emp.n_p_context_atx(find(outp_atx.p_context_atx_p<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_context_atx_n<alpha1),emp.n_n_context_atx(find(outp_atx.p_context_atx_n<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_context_atx_p<alpha2),emp.n_p_context_atx(find(outp_atx.p_context_atx_p<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_context_atx_n<alpha2),emp.n_n_context_atx(find(outp_atx.p_context_atx_n<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_context_atx_p<alpha3),emp.n_p_context_atx(find(outp_atx.p_context_atx_p<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_context_atx_n<alpha3),emp.n_n_context_atx(find(outp_atx.p_context_atx_n<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
tp_editplots
pos(5,:)=get(gca,'Position')
axis([freq_start-1 nfreq -0.25 .25])

subplot(4,3,8); hold on
plot(freq_start:nfreq,emp.n_n_context_dpz(freq_start:nfreq),'b-','linewidth',2)
plot(freq_start:nfreq,emp.n_p_context_dpz(freq_start:nfreq),'r-','linewidth',2)
axis([freq_start-1 nfreq lims_ctx(2,1) lims_ctx(2,end)])
set(gca,'tickdir','out','ytick',[-0.25 0 0.25],'yticklabel',num2cell([-25 0 25]))
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
xlabel('Carrier frequency [Hz]'); 
plot(find(outp_atx.p_context_dpz_p<alpha1),emp.n_p_context_dpz(find(outp_atx.p_context_dpz_p<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_context_dpz_n<alpha1),emp.n_n_context_dpz(find(outp_atx.p_context_dpz_n<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.p_context_dpz_p<alpha2),emp.n_p_context_dpz(find(outp_atx.p_context_dpz_p<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_context_dpz_n<alpha2),emp.n_n_context_dpz(find(outp_atx.p_context_dpz_n<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.p_context_dpz_p<alpha3),emp.n_p_context_dpz(find(outp_atx.p_context_dpz_p<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.p_context_dpz_n<alpha3),emp.n_n_context_dpz(find(outp_atx.p_context_dpz_n<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
tp_editplots
pos(6,:)=get(gca,'Position')
axis([freq_start-1 nfreq -0.25 .25])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_lineplots_allfreqs_stats_v%d.pdf',v));

%% (5) PLOT ALTERED CORRELATIONS (PER VOXEL)

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

para          = [];
para.nfreq    = [1:size(emp.n_p_atx,1)]; % freqs 1 - 13
para.alpha    = 0.05; % alpha for adjacency
para.ver      = v;
para.nperm    = 10000;
para.nsubs    = 200;
para.type     = 'local';
para.allperms = para.nperm/para.nsubs;
para.emp      = emp;
para.cleaned  = 0;
para.correction_method = 'single_threshold';
% ---------
% Obtain stats for atomoxetine condition

outp    = pupmod_src_powcorr_getstatistics(para);
% ---------
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
%% PLOT RESULTS FOR ATX (TASK) AND DPZ (REST)

close all
for ifoi = 6:9
  icond = 2;

  cmap  = autumn;
  par   = nanmean(emp.n_p_atx_pervoxel(:,ifoi,icond),2);

  par(mean(outp.pval_vox_p_atx(:,icond,ifoi),3)>=0.05) = 0;

  cmap      = [cmap; 0.98*ones(1,3); cmap];
  para      = [];
  para.clim = [-0.4 0.4];
  para.cmap = cmap;
  para.grid = grid;
  para.dd   = 0.75;
  para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_atx_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
  tp_plot_surface(par,para)
  
end
%% MAP DONEPEZIL EFFECT

for ifoi =  6:10
  
  icond     = 1;
  cmap      = autumn;
  cmap(:,1) = 0; 
  cmap(:,3) = 1;

  par = nanmean(emp.n_n_dpz_pervoxel(:,ifoi,icond),2);

  par(outp.pval_vox_n_dpz(:,icond,ifoi)>=0.05) = 0;
% par(par<prctile(par,95))=0
  cmap      = [cmap(:,:); 0.98*ones(1,3); cmap(:,:)];
  para      = [];
  para.clim = [-0.4 0.4];
  para.cmap = cmap;
  para.grid = grid;
  para.dd   = 0.75;
  para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_dpz_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
  tp_plot_surface(par,para);
end


%% LOAD SECOND LEVEL STATS (POOLED ACROSS SIGNIFICANT FREQS)
% -----------------------------------------
% load second level permutation test results
% this accounts for double dipping (see hawellek et al. 2014)

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

% ---------
% Obtain stats for atomoxetine condition
para.cond     = 'atx';
% fc(:,:,:,:,:,18) = nanmean(fc(:,:,:,:,:,6:9),6);
% dpz effect
% fc(:,:,:,:,:,19) =  nanmean(fc(:,:,:,:,:,6:10),6);
% fc(:,:,:,:,:,20) =  nanmean(fc(:,:,:,:,:,4:6),6);

para          = [];
para.alpha    = 0.05; % alpha for adjacency
para.ver      = v;
para.nperm    = 10000;
para.nsubs    = 200;
para.type     = 'global';
para.cond     = 'atx';
para.allperms = para.nperm/para.nsubs;
para.nfreq    = 1:17;
para.emp      = pupmod_compute_altered_correlations(fc,para);
para.cleaned  = 0;
para.correction_method = 'ranks';
para.v = v;

[outp_atx]    = pupmod_src_powcorr_getstatistics(para);

emp=para.emp;
% % ---------
%% ATOMOXETINE
% --------
ifoi =1  ; icond = 2;

par = nanmean(emp.n_p_atx_pervoxel(:,ifoi,icond),2);
par(outp.pval_vox_p_atx(:,icond,ifoi)>=0.05) = 0;

cmap = autumn;

cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_2nd_atx_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para)

%% DONEPEZIL
% --------
ifoi =2; icond = 1;

par = nanmean(emp.n_n_dpz_pervoxel(:,ifoi,icond),2);
par(outp.pval_vox_n_dpz(:,icond,ifoi)>=0.05) = 0;

cmap = autumn; cmap(:,1) = 0; cmap(:,3) = 1;

cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_2nd_dpz_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para)

%% TASK VS REST MAP

ifoi =3; icond = 1;

par = nanmean(emp.taskvsrest_n_pervoxel(:,ifoi),2);
par(outp.pval_vox_n_tvr(:,icond,ifoi)>=0.01) = 0;

cmap = autumn; cmap(:,1) = 0; cmap(:,3) = 1;

cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_2nd_tvr_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para)



%% TASK VS REST
% TASK VS REST: lineplots
% ----------------------

% Obtain stats for task vs rest
para.cond     = 'taskvsrest';
[outp_tvr]    = pupmod_src_powcorr_getstatistics(para);

%%
figure; set(gcf,'color','w');
lims = [0 0.35];
lims_lab = num2cell([0 35]);
lims_ctx = [-0.25 0 0.30; -0.25 0 0.25];

markersize = 4;
  
alpha1 = 0.05;
alpha2 = 0.01;
alpha3 = 0.001;
subplot(4,2,2); hold on
plot(emp.taskvsrest_p,'r-','linewidth',2)
plot(emp.taskvsrest_n,'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',num2cell([4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
axis([0 17 lims(1)-0.05 lims(end)])
plot(find(outp_atx.pval_taskvsrest_p_corr<alpha1),emp.taskvsrest_p(find(outp_atx.pval_taskvsrest_p_corr<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.pval_taskvsrest_n_corr<alpha1),emp.taskvsrest_n(find(outp_atx.pval_taskvsrest_n_corr<alpha1)),'ko','markersize',markersize,'markerfacecolor','k')
plot(find(outp_atx.pval_taskvsrest_p_corr<alpha2),emp.taskvsrest_p(find(outp_atx.pval_taskvsrest_p_corr<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.pval_taskvsrest_n_corr<alpha2),emp.taskvsrest_n(find(outp_atx.pval_taskvsrest_n_corr<alpha2)),'ko','markersize',markersize,'markerfacecolor','w')
plot(find(outp_atx.pval_taskvsrest_p_corr<alpha3),emp.taskvsrest_p(find(outp_atx.pval_taskvsrest_p_corr<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
plot(find(outp_atx.pval_taskvsrest_n_corr<alpha3),emp.taskvsrest_n(find(outp_atx.pval_taskvsrest_n_corr<alpha3)),'ko','markersize',markersize,'markerfacecolor','m')
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
tp_editplots
pos(2,:)=get(gca,'Position');

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_tvr_lineplots_allfreqs_stats_v%d.pdf',v));

%%
ifoi = 4:6; icond = 1;

cmap = autumn;
cmap(:,1) = 0; cmap(:,3) = 1;

par = nanmean(emp.taskvsrest_n_pervoxel(:,ifoi),2);
% par(outp_tvr.pval_p_tvr(:,ifoi(1))>=0.05) = 0;
% par(par<prctile(par,95))=0
cmap = [cmap(:,:); 0.98*ones(1,3); cmap(:,:)];
para      = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_tvr_nomask_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para)


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
axis([0 size(emp.n_p_atx,1) lims(1)-0.05 lims(end)+0.175])
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

%% COMPUTE ALTERED CORRELATIONS FOR DIFFERENT ALPHA LEVELS

clear emp;

alphas = [0.01:0.01:0.10]
for ialpha = 1 : 10
  
  para.alpha = alphas(ialpha);
  para.nfreq = 1 : 17;
  
  emp{ialpha} = pupmod_compute_altered_correlations(fc,para);
  
end

save(sprintf('~/pupmod/proc/conn/emp_severalalphas_v%d.mat',v),'emp')

%% PLOT ALTERED CORRELATION FOR DIFFERENT ALPHA LEVELS
v = 3;
load ~/pupmod/proc/conn/emp_severalalphas_v3.mat

figure; set(gcf,'color','w');
cmap = cbrewer('seq', 'Greys', 15,'pchip');
cmap = cmap(4:13,:);

subplot(4,3,1);hold on
for i = 1 : 10
  plot(100*emp{i}.n_p_atx(:,1),'color',cmap(i,:))
end
axis([0 17 -5 40])

set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')
subplot(4,3,4);hold on
for i = 1 : 10
  plot(100*emp{i}.n_p_atx(:,2),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')

subplot(4,3,2);hold on
for i = 1 : 10
  plot(100*emp{i}.n_n_atx(:,1),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')
subplot(4,3,5);hold on
for i = 1 : 10
  plot(100*emp{i}.n_n_atx(:,2),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]'); %ylabel('Fraction of significantly altered correlations [%]')

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_alteredcorr_alternative_alphas_atx_v%d.eps',v))

figure; set(gcf,'color','w');
cmap = cbrewer('seq', 'Greys', 15,'pchip');
cmap = cmap(4:13,:);

subplot(4,3,1);hold on
for i = 1 : 10
  plot(100*emp{i}.n_p_dpz(:,1),'color',cmap(i,:))
end
axis([0 17 -5 40])

set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')
subplot(4,3,4);hold on
for i = 1 : 10
  plot(100*emp{i}.n_p_dpz(:,2),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')

subplot(4,3,2);hold on
for i = 1 : 10
  plot(100*emp{i}.n_n_dpz(:,1),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]')
subplot(4,3,5);hold on
for i = 1 : 10
  plot(100*emp{i}.n_n_dpz(:,2),'color',cmap(i,:))
end
axis([0 17 -5 40])
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
tp_editplots; xlabel('Carrier frequency [Hz]'); %ylabel('Fraction of significantly altered correlations [%]')

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_alteredcorr_alternative_alphas_dpz_v%d.eps',v))


%% PLOT FC IN YEO ATLAS SPACE
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
grid = sa.grid_cortex_lowres;

lab=tp_aal2yeo(grid);
cmap = cbrewer('qual', 'Dark2', 7); 
para = [];
para.nfreq = 1:17;
para.alpha = 0.05;

for i = 1 : 7
emp{i} = pupmod_compute_altered_correlations(fc(lab==i,lab==i,:,:,:,:),para);
end

%%

para          = [];
para.nfreq    = 1:17; % freqs 1 - 13
para.alpha    = 0.05; % alpha for adjacency
para.ver      = v;
para.nperm    = 10000;
para.nsubs    = 200;
para.type     = 'global';
para.allperms = para.nperm/para.nsubs;
para.emp      = emp;
para.cleaned  = 0;
para.correction_method = 'single_threshold';
% ---------
% Obtain stats for atomoxetine condition
para.cond     = 'atx';
[outp_atx]    = pupmod_src_powcorr_getstatistics_yeo(para);


%%
cmap = cbrewer('qual', 'Dark2', 7);
figure; set(gcf,'color','w'); hold on
subplot(4,2,1); hold on
alpha1 = 0.05;
markersize=3;
for i = [1 2 3 4 6 7]

  plot(100*emp{i}.n_p_atx(1:17,2),'color',cmap(i,:))
  [cl,n] = bwlabeln(outp_atx(i).p_cnt1_p<alpha1);
  
  plot(find(cl),ones(1,sum(cl)).*100.*emp{i}.n_p_atx(find(cl),2),'.','markersize',12,'color',cmap(i,:))
%  
%   for in = 1 : n
%     if sum(cl==in)>1
%       line([min(find(cl==in)) max(find(cl==in))],[35+2*i 35+2*i],'color',cmap(i,:))
%     else
%       plot(find(cl==n),35+2*i,'ko','markersize',markersize,'markerfacecolor',cmap(i,:),'markeredgecolor',cmap(i,:))
%     end
%   end
end
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]')
tp_editplots
axis([0 17 0 40])
subplot(4,2,2); hold on
for i = [1 2 3 4 6 7]

   plot(100*emp{i}.n_n_dpz(1:17,1),'color',cmap(i,:))
  [cl,n] = bwlabeln(outp_atx(i).p_res2_n<alpha1);
 
    plot(find(cl),ones(1,sum(cl)).*100.*emp{i}.n_n_dpz(find(cl),1),'.','markersize',12,'color',cmap(i,:))

  
end
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]')
tp_editplots
axis([0 17 0 40])

subplot(4,2,3); hold on
for i = [1 2 3 4 6 7]

  plot(100*emp{i}.n_n_dpz(1:17,1),'color',cmap(i,:))
end

legend('Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','DMN')

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_yeoatlas_within_v%d.eps',v))
%% PLOT YEO ATLAS
cols = {'Greens';'Oranges';'Blues';'RdPu';'YlGn';'Reds'};
for i = 1:7
  if i == 7
    cmap = cbrewer('div','BrBG',256); cmap=cmap(1:128,:); cmap=cmap(end:-1:1,:);
  else
    cmap = cbrewer('seq',cols{i},128);

  end
    cmap = [cmap; 0.98*ones(1,3); cmap];

par = lab==i;

% par(outp_atx.pval_p_atx(:,icond,ifoi)>=0.05) = 0;
% par(par<0.4)=0;

para = [];
para.clim = [-0.1 0.1];
para.cmap = cmap;
para.grid = grid;
para.dd = 1;
para.fn = sprintf('~/pupmod/plots/yeo_reg%d.png',i);
tp_plot_surface(par,para)

end

%% BETWEEN YEO REGIONS
clear emp;
para = [];
para.nfreq = 1:25;
para.alpha = 0.01;
for i = 1 : 7
  for j = 1 : 7
    idx1=find(lab==i);
    idx2=find(lab==j);
    
     fc_tmp = fc(idx1,idx2,:,:,:,:);
      
     emp{i}{j} = pupmod_compute_altered_correlations(fc_tmp,para);

  end
end

within = zeros(25,1); across = zeros(25,1); w= 0; a=0;
for i = 1 : 7
  for j = 1 : 7
    if i == j
      w = w + 1;
      within = within + 100*emp{i}{j}.n_n_dpz(:,1);
    else
      a = a + 1;
      across = across + 100*emp{i}{j}.n_n_dpz(:,1);
    end
  end
end

within=within./w;
across=across./a;
%% PLOT MATRIX FOR SPECIFIED FREQ

foi_atx = [10:12];
foi_dpz = [13 14];

for i = 1 : 7
  for j = 1 : 7
    kk_atx(i,j) = mean(100*emp{i}{j}.n_p_atx(foi_atx,2),1);
    kk_dpz(i,j) = mean(100*emp{i}{j}.n_n_dpz(foi_dpz,1),1);

  end
end
    
figure; set(gcf,'color','w');

subplot(1,2,1);
imagesc(kk_atx,[0 30])
colormap(plasma)
axis square; tp_editplots
set(gca,'tickdir','out','xtick',[1:7],'xticklabel',{'Visual';'Somatomotor';'Dorsal Att.';'Ventral Att.';'Limbic';'Frontoparietal';'DMN'})
xtickangle(45)
set(gca,'tickdir','out','ytick',[1:7],'yticklabel',{'Visual';'Somatomotor';'Dorsal Att.';'Ventral Att.';'Limbic';'Frontoparietal';'DMN'})
tp_colorbar()

subplot(1,2,2);
imagesc(kk_dpz,[0 30])
colormap(plasma)
axis square; tp_editplots
set(gca,'tickdir','out','xtick',[1:7],'xticklabel',{'Visual';'Somatomotor';'Dorsal Att.';'Ventral Att.';'Limbic';'Frontoparietal';'DMN'})
xtickangle(45)
set(gca,'tickdir','out','ytick',[1:7],'yticklabel',{'Visual';'Somatomotor';'Dorsal Att.';'Ventral Att.';'Limbic';'Frontoparietal';'DMN'})
tp_colorbar()

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_yeoatlas_across_fatx%s_fdpz%s_v%d.eps',regexprep(num2str(foi_atx),' ',''),regexprep(num2str(foi_dpz),' ',''),v))


