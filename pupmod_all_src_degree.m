%% pupmod_all_src_degree
% plot degree of cleaned signal
% obtain cleanined signal from pupmod_all_src_peripheral*** *check)
% Goal: replicate hipp et al., nn and plot stats



clear
% version: 12 coarse cortex, 1 AAL

v = 12;
outdir = '~/pupmod/proc/conn/';

cleandat = pupmod_loadpowcorr(v);
%% COMPUTE "RELATIVE" DEGREE
% --------------------------------

fcsize = size(cleandat,1);
para = [];
para.alpha = 0.01;
para.nfreq = 13;

% Atomoxetine (icond = 2) during Rest 
deg_atx = tp_degree(cleandat(:,:,:,[1 2],1,:),para);
deg_atx_vox = squeeze(nansum(deg_atx)/fcsize);
deg_atx = squeeze(nansum(reshape(deg_atx,[fcsize^2 13 2]))/fcsize^2);
% Atomoxetine (icond = 3) during Task 
deg_atx_task = tp_degree(cleandat(:,:,:,[1 2],2,:),para);
deg_atx_task_vox = squeeze(nansum(deg_atx_task)/fcsize);
deg_atx_task = squeeze(nansum(reshape(deg_atx_task,[fcsize^2 13 2]))/fcsize^2);
% ------------------
% Donepezil (icond = 3) during Rest 
deg_dpz = tp_degree(cleandat(:,:,:,[1 3],1,:),para);
deg_dpz_vox = squeeze(nansum(deg_dpz)/fcsize);
deg_dpz = squeeze(nansum(reshape(deg_dpz,[fcsize^2 13 2]))/fcsize^2);
% Donepezil (icond = 3) during Task 
deg_dpz_task = tp_degree(cleandat(:,:,:,[1 3],2,:),para);
deg_dpz_task_vox = squeeze(nansum(deg_dpz_task)/fcsize);
deg_dpz_task = squeeze(nansum(reshape(deg_dpz_task,[fcsize^2 13 2]))/fcsize^2);
% ------------------

%% PLOT  RESULTS W/O STATISTICS
% --------------------------------

foi_range = unique(round(2.^[1:0.5:7]));

figure; set(gcf,'color','w'); hold on
subplot(2,2,1); hold on

area(deg_atx(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
k = area(deg_atx(:,2),'facecolor',[1 0.5 0.2],'edgecolor',[0.8 0.8 0.8]);
alpha(k,0.6)
set(gca,'xtick',1:2:13,'xticklabel',foi_range(1:2:13))
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([1 13 0 0.5]);

subplot(2,2,2); hold on

area(deg_atx_task(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
k = area(deg_atx_task(:,2),'facecolor',[1 0.5 0.2],'edgecolor',[0.8 0.8 0.8]);
alpha(k,0.6)
set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([1 13 0 0.5]);

subplot(2,2,3); hold on

area(deg_dpz(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
k = area(deg_dpz(:,2),'facecolor',[0.2 0.5 1],'edgecolor',[0.8 0.8 0.8]);
alpha(k,0.6)
set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([1 13 0 0.5]);

subplot(2,2,4); hold on

area(deg_dpz_task(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
k = area(deg_dpz_task(:,2),'facecolor',[0.2 0.5 1],'edgecolor',[0.8 0.8 0.8]);
alpha(k,0.6)
set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([1 13 0 0.5]);

%% CONCATENATE PERMUTATONS AND COMPUTE STATS
v = 12;
nperm = 50000; 
par.subs = 100;
par.allperms = nperm/par.subs;

for iperm = 1 : par.allperms
  iperm
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  
  perm_k_atx(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_atx;
  perm_k_atx_vox(:,:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_atx_pervoxel;
  perm_k_dpz(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs)   = outp.perm_k_dpz;
  perm_k_dpz_vox(:,:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_dpz_pervoxel;

end

clear outp

for ifoi = 1 : 13
  
  ifoi 
  p_atx_task(ifoi) = 1-sum(abs(deg_atx_task(ifoi,2)-deg_atx_task(ifoi,1)) > abs(squeeze(perm_k_atx(ifoi,2,2,:)-perm_k_atx(ifoi,1,2,:))))/50000;
  p_atx_rest(ifoi) = 1-sum(abs(deg_atx(ifoi,2)-deg_atx(ifoi,1)) > abs(squeeze(perm_k_atx(ifoi,2,1,:)-perm_k_atx(ifoi,1,1,:))))/50000;
  p_dpz_rest(ifoi) = 1-sum(abs(deg_dpz(ifoi,2)-deg_dpz(ifoi,1)) > abs(squeeze(perm_k_dpz(ifoi,2,1,:)-perm_k_dpz(ifoi,1,1,:))))/50000;
  p_dpz_task(ifoi) = 1-sum(abs(deg_dpz_task(ifoi,2)-deg_dpz_task(ifoi,1)) > abs(squeeze(perm_k_dpz(ifoi,2,2,:)-perm_k_dpz(ifoi,1,2,:))))/50000;
  
  p_atx_vox_task(:,ifoi) = 1-sum(abs(deg_atx_task_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1)) > abs(squeeze(perm_k_atx_vox(:,ifoi,2,2,:)-perm_k_atx_vox(:,ifoi,1,2,:))),2)/50000;
  p_atx_vox_rest(:,ifoi) = 1-sum(abs(deg_atx_vox(:,ifoi,2)-deg_atx_vox(:,ifoi,1)) > abs(squeeze(perm_k_atx_vox(:,ifoi,2,1,:)-perm_k_atx_vox(:,ifoi,1,1,:))),2)/50000;
  p_dpz_vox_rest(:,ifoi) = 1-sum(abs(deg_dpz_vox(:,ifoi,2)-deg_dpz_vox(:,ifoi,1)) > abs(squeeze(perm_k_dpz_vox(:,ifoi,2,1,:)-perm_k_dpz_vox(:,ifoi,1,1,:))),2)/50000;
  p_dpz_vox_task(:,ifoi) = 1-sum(abs(deg_dpz_task_vox(:,ifoi,2)-deg_dpz_task_vox(:,ifoi,1)) > abs(squeeze(perm_k_dpz_vox(:,ifoi,2,2,:)-perm_k_dpz_vox(:,ifoi,1,2,:))),2)/50000;
  

end

clear perm_k_atx perm_k_dpz perm_k_atx_vox perm_k_dpz_vox
%% PLOT ON SURFACE
v = 12
ifoi = 13;
% icond = 2;

% load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],2,6,v));

var2plot_atx(:,1) = (deg_atx_vox(:,ifoi,2)-deg_atx_vox(:,ifoi,1)).*(p_atx_vox_rest(:,ifoi)<0.01);
var2plot_atx(:,2) = (deg_atx_task_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1)).*(p_atx_vox_task(:,ifoi)<0.01);
var2plot_dpz(:,1) = (deg_dpz_vox(:,ifoi,2)-deg_dpz_vox(:,ifoi,1)).*(p_dpz_vox_rest(:,ifoi)<0.01);
var2plot_dpz(:,2) = (deg_dpz_task_vox(:,ifoi,2)-deg_dpz_task_vox(:,ifoi,1)).*(p_dpz_vox_task(:,ifoi)<0.05);

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

vc    = sa_template.vc;
g1    = grid;
g2    = sa_template.cortex10K.vc;
dd    = .5;

interp_atx(:,1) = spatfiltergauss(var2plot_atx(:,1),g1,dd,g2);
interp_atx(:,2) = spatfiltergauss(var2plot_atx(:,2),g1,dd,g2);
interp_dpz(:,1) = spatfiltergauss(var2plot_dpz(:,1),g1,dd,g2);
interp_dpz(:,2) = spatfiltergauss(var2plot_dpz(:,2),g1,dd,g2);

%% PLOT DEGREE, PLACEBO
% --------------------------------

ifoi = 7; icond = 1;

cmap = plasma;

par = deg_atx_vox(:,ifoi,1);

cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
tp_plot_surface(par,sa_template,para)
