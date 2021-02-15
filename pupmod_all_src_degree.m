%% pupmod_all_src_degree
% plot degree of cleaned signal
% obtain cleanined signal from pupmod_all_src_peripheral*** *check)
% Goal: replicate hipp et al., nn and plot stats



clear
% version: 12 coarse cortex, 1 AAL

v = 3;
outdir = '~/pupmod/proc/conn/';

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
cleandat = pupmod_loadpowcorr(v,SUBJLIST,0);
%% COMPUTE "RELATIVE" DEGREE
% --------------------------------
fcsize = size(cleandat,1);
mask    = logical(triu(ones(400,400),1));
para = [];
para.alpha = 0.01;
para.nfreq = 17;
para.absolute = 0;
para.relative_degree = 0;
para.clustering = 1;
para.transitivity = 1;

% Atomoxetine (icond = 2) during Rest 
outp_atx_rest = tp_degree(cleandat(:,:,:,[1 2],1,:),para);
% Atomoxetine (icond = 2) during Task 
outp_atx_task = tp_degree(cleandat(:,:,:,[1 2],2,:),para);
% ------------------
% Donepezil (icond = 3) during Rest 
outp_dpz_rest = tp_degree(cleandat(:,:,:,[1 3],1,:),para);
% Donepezil (icond = 3) during Task 
outp_dpz_task = tp_degree(cleandat(:,:,:,[1 3],2,:),para);
% ------------------

%% PLOT  RESULTS W/O STATISTICS
% --------------------------------

foi_range       = 2.^[1:.25:6];
figure; set(gcf,'color','w'); hold on
subplot(3,2,1); hold on
title('Atx - Rest')

% area(deg_atx(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
plot(outp_atx_rest.tot_degree(:,1),'k')
% plot(outp_atx_task.tot_degree(:,1),'k','linestyle',':');
% area(deg_atx_task(:,1),'facecolor',[0.9 0.9 0.9],'edgecolor',[1 1 1]);
% plot(deg_atx_task(:,1),'color',[0 0 0],'linestyle',':')
% k = area(deg_atx(:,2),'facecolor',[1 0.5 0.2],'edgecolor',[0.8 0.8 0.8]);
% alpha(k,0.6)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([0 length(foi_range) 0 30]);

subplot(3,2,2); hold on
title('Atx - Task')
plot(outp_atx_task.tot_degree(:,1),'k')
plot(outp_atx_task.tot_degree(:,2),'k','linestyle',':');
% k = area(deg_atx_task(:,2),'facecolor',[1 0.5 0.2],'edgecolor',[0.8 0.8 0.8]);
% alpha(k,0.6)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([0 length(foi_range) 0 30]);

subplot(3,2,3); hold on
title('Dpz - Rest')
plot(outp_dpz_rest.tot_degree(:,1),'k')
plot(outp_dpz_rest.tot_degree(:,2),'k','linestyle',':');
% area(deg_dpz(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
% k = area(deg_dpz(:,2),'facecolor',[0.2 0.5 1],'edgecolor',[0.8 0.8 0.8]);
% alpha(k,0.6)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([0 length(foi_range) 0 40]);

subplot(3,2,4); hold on
title('Dpz - Task')
plot(outp_dpz_task.tot_degree(:,1),'k')
plot(outp_dpz_task.tot_degree(:,2),'k','linestyle',':');
% area(deg_dpz_task(:,1),'facecolor',[0.8 0.8 0.8],'edgecolor',[1 1 1]);
% k = area(deg_dpz_task(:,2),'fvacecolor',[0.2 0.5 1],'edgecolor',[0.8 0.8 0.8]);
% alpha(k,0.6)
set(gca,'tickdir','out','xtick',[1 5 9 13 17],'xticklabel',[4 8 16 32 64])
xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
tp_editplots
axis([0 length(foi_range) 0 40]);

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_degree_abs%d_rel%d_v%d.pdf',para.absolute,para.relative_degree,v))
%% CONCATENATE PERMUTATONS AND COMPUTE STATS
v = 23;
nperm = 10000; 
par.subs = 250;
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

for ifoi = 1 : 25
  
  ifoi 
  p_atx_task(ifoi) = 1-sum(abs(deg_atx_task(ifoi,2)-deg_atx_task(ifoi,1)) > abs(squeeze(perm_k_atx(ifoi,2,2,:)-perm_k_atx(ifoi,1,2,:))))/nperm;
  p_atx_rest(ifoi) = 1-sum(abs(deg_atx(ifoi,2)-deg_atx(ifoi,1)) > abs(squeeze(perm_k_atx(ifoi,2,1,:)-perm_k_atx(ifoi,1,1,:))))/nperm;
  p_dpz_rest(ifoi) = 1-sum(abs(deg_dpz(ifoi,2)-deg_dpz(ifoi,1)) > abs(squeeze(perm_k_dpz(ifoi,2,1,:)-perm_k_dpz(ifoi,1,1,:))))/nperm;
  p_dpz_task(ifoi) = 1-sum(abs(deg_dpz_task(ifoi,2)-deg_dpz_task(ifoi,1)) > abs(squeeze(perm_k_dpz(ifoi,2,2,:)-perm_k_dpz(ifoi,1,2,:))))/nperm;
  
  p_atx_vox_task(:,ifoi) = 1-sum(abs(deg_atx_task_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1)) > abs(squeeze(perm_k_atx_vox(:,ifoi,2,2,:)-perm_k_atx_vox(:,ifoi,1,2,:))),2)/nperm;
  p_atx_vox_rest(:,ifoi) = 1-sum(abs(deg_atx_vox(:,ifoi,2)-deg_atx_vox(:,ifoi,1)) > abs(squeeze(perm_k_atx_vox(:,ifoi,2,1,:)-perm_k_atx_vox(:,ifoi,1,1,:))),2)/nperm;
  p_dpz_vox_rest(:,ifoi) = 1-sum(abs(deg_dpz_vox(:,ifoi,2)-deg_dpz_vox(:,ifoi,1)) > abs(squeeze(perm_k_dpz_vox(:,ifoi,2,1,:)-perm_k_dpz_vox(:,ifoi,1,1,:))),2)/nperm;
  p_dpz_vox_task(:,ifoi) = 1-sum(abs(deg_dpz_task_vox(:,ifoi,2)-deg_dpz_task_vox(:,ifoi,1)) > abs(squeeze(perm_k_dpz_vox(:,ifoi,2,2,:)-perm_k_dpz_vox(:,ifoi,1,2,:))),2)/nperm;
  

end

% clear perm_k_atx perm_k_dpz perm_k_atx_vox perm_k_dpz_vox
%% PLOT ON SURFACE
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
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
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

for ifoi = 10:12
  for icond = 1 :2
% ifoi = 1; icond = 1;

cmap = plasma;
% cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
par = deg_atx_vox(:,ifoi,icond);

% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);

% cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [0.1 0.5]
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_degree_atx_rest_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para);

par = deg_atx_task_vox(:,ifoi,icond);

% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);

% cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [0.1 0.5]
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_degree_atx_task_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para);


par = deg_atx_task_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [-0.3 0.3]
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_degree_atx_task_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para);


par = deg_dpz_vox(:,ifoi,icond);

% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);

% cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [0.1 0.5]
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_degree_dpz_rest_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para);

par = deg_dpz_task_vox(:,ifoi,icond);

% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
% par = deg_atx_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);

% cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [0.1 0.5]
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn = sprintf('~/pupmod/plots/pupmod_plot_degree_dpz_task_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
tp_plot_surface(par,para);





  end
  close all
end
par = deg_atx_vox(:,ifoi,2) - deg_atx_vox(:,ifoi,2);

