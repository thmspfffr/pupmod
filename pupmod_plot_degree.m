%% pupmod_all_src_degree
% --------
% plot degree of cleaned signal
% obtain cleanined signal from pupmod_all_src_peripheral***
% Goal: replicate hipp et al., nn and plot stats
% --------
% Permutation distribution computed in pupmod_all_src_degree_permtest.m



clear
% version: 12 coarse cortex, 1 AAL
addpath ~/pupmod/matlab
v = 3;
outdir = '~/pupmod/proc/conn/';
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
fc = pupmod_loadpowcorr(v,SUBJLIST,1);

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
%% COMPUTE "RELATIVE" DEGREE
% --------------------------------

fcsize               = size(fc,1);
para                 = [];
para.alpha           = 0.01;
para.nfreq           = 17;
para.absolute        = 0;
para.relative_degree = 0;

% Atomoxetine (icond = 2) during Rest 
deg_atx = tp_degree(fc(:,:,:,[1 2],1,:),para);
% deg_atx_vox = squeeze(nansum(deg_atx.node_degree)/fcsize);
% deg_atx = squeeze(nansum(reshape(deg_atx,[fcsize^2 para.nfreq 2]))/fcsize^2);
% Atomoxetine (icond = 3) during Task 
deg_atx_task = tp_degree(fc(:,:,:,[1 2],2,:),para);
% deg_atx_task_vox = squeeze(nansum(deg_atx_task)/fcsize);
% deg_atx_task = squeeze(nansum(reshape(deg_atx_task,[fcsize^2 para.nfreq 2]))/fcsize^2);
% ------------------
% Donepezil (icond = 3) during Rest 
deg_dpz = tp_degree(fc(:,:,:,[1 3],1,:),para);
% deg_dpz_vox = squeeze(nansum(deg_dpz)/fcsize);
% deg_dpz = squeeze(nansum(reshape(deg_dpz,[fcsize^2 para.nfreq 2]))/fcsize^2);
% Donepezil (icond = 3) during Task 
deg_dpz_task = tp_degree(fc(:,:,:,[1 3],2,:),para);
% deg_dpz_task_vox = squeeze(nansum(deg_dpz_task)/fcsize);
% deg_dpz_task = squeeze(nansum(reshape(deg_dpz_task,[fcsize^2 para.nfreq 2]))/fcsize^2);
% ------------------

%% PLOT  RESULTS W/O STATISTICS
% --------------------------------

foi_range = 2.^(2:.25:6);
fc_avg = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1:2,1,:),1),2),3),4));

figure; set(gcf,'color','w'); hold on
subplot(4,2,1); hold on
yyaxis left
plot(deg_atx.tot_degree(:,1),'color',[0 0 0],'linewidth',2);
yyaxis right 
plot(fc_avg,'color',[0.8 0.8 0.8],'linewidth',2);

% plot(deg_atx.tot_degree(:,2),'color',[1 0.5 0.2],'linewidth',2);
yyaxis left
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
ylabel('Degree [%]')
axis(gca,[0 18 -5 30]);

yyaxis right
set(gca,'ytick',[0 0.05 0.1 0.15],'yticklabel',[0 0.05 0.10 0.15])
axis(gca,[0 18 -0.01 0.06]);

xlabel('Carrier frequency [Hz]'); ylabel('Mean FC')
tp_editplots

line([9 9],[0 30])
line([6 6],[0 0.15])

% subplot(4,2,2); hold on
% 
% plot(deg_atx_task.tot_degree(:,1),'color',[0.8 0.8 0.8],'linewidth',2);
% plot(deg_atx_task.tot_degree(:,2),'color',[1 0.5 0.2],'linewidth',2);
% set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
% xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
% tp_editplots
% axis([0 14 -0.02 0.6]);
% 
% subplot(4,2,3); hold on
% 
% plot(deg_dpz.tot_degree(:,1),'color',[0.8 0.8 0.8],'linewidth',2);
% plot(deg_dpz.tot_degree(:,2),'color',[0.2 0.5 1],'linewidth',2);
% set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
% xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
% tp_editplots
% axis([0 14 -0.02 0.6]);
% 
% subplot(4,2,4); hold on
% 
% plot(deg_dpz_tas.tot_degreek(:,1),'color',[0.8 0.8 0.8],'linewidth',2);
% plot(deg_dpz_task.tot_degree(:,2),'color',[0.2 0.5 1],'linewidth',2);
% set(gca,'xtick',(1:2:13),'xticklabel',foi_range(1:2:13))
% xlabel('Carrier frequency [Hz]'); ylabel('Degree [%]')
% tp_editplots
% axis([0 14 -0.02 0.6]);

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_degree_lineplots_rel%d_v%d.pdf',0,v));

%% LOAD AND GENERATE PERMUTATION DISTRIBUTION
% Details defined in par structure

v = 3;
clear atx dpz outp

par           = [];
nperm         = 10000; 
par.subs      = 250;
par.allperms  = nperm/par.subs;

for iperm = 1 : par.allperms
  fprintf('Loading permutation file %d/%d (ATX) ...\n',iperm,par.allperms)
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  atx.perm_k_atx(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs)     = single(outp.perm_k_atx);
  atx.perm_k_atx_vox(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = squeeze(outp.perm_k_atx_pervoxel(:,:,2,:,:)-outp.perm_k_atx_pervoxel(:,:,1,:,:));
  clear outp
end

for iperm = 1 : par.allperms
  fprintf('Loading permutation file %d/%d (DPZ) ...\n',iperm,par.allperms)
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  dpz.perm_k_dpz(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs)     = single(outp.perm_k_dpz);
  dpz.perm_k_dpz_vox(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = squeeze(outp.perm_k_dpz_pervoxel(:,:,2,:,:)-outp.perm_k_dpz_pervoxel(:,:,1,:,:));
  clear outp
end

%% Obtain p-values (corrected using single threshold permutation - across space)
% See Nichols & Holmes (2002) Hum Brain Mapp for details
clear p* max*

for ifoi = 1 : 13
  
  % Compute global p-values
  p_atx_task(ifoi) = 1-sum(abs(deg_atx_task(ifoi,2)-deg_atx_task(ifoi,1)) > abs(squeeze(atx.perm_k_atx(ifoi,2,2,:)-atx.perm_k_atx(ifoi,1,2,:))))/nperm;
  p_atx_rest(ifoi) = 1-sum(abs(deg_atx(ifoi,2)-deg_atx(ifoi,1)) > abs(squeeze(atx.perm_k_atx(ifoi,2,1,:)-atx.perm_k_atx(ifoi,1,1,:))))/nperm;
  p_dpz_rest(ifoi) = 1-sum(abs(deg_dpz(ifoi,2)-deg_dpz(ifoi,1)) > abs(squeeze(dpz.perm_k_dpz(ifoi,2,1,:)-dpz.perm_k_dpz(ifoi,1,1,:))))/nperm;
  p_dpz_task(ifoi) = 1-sum(abs(deg_dpz_task(ifoi,2)-deg_dpz_task(ifoi,1)) > abs(squeeze(dpz.perm_k_dpz(ifoi,2,2,:)-dpz.perm_k_dpz(ifoi,1,2,:))))/nperm;
  
  % Compute p-values for each voxel
  if 0
    % This is already the difference, computed one section earlier  
    max_atx(:,1) = squeeze(max(squeeze(atx.perm_k_atx_vox(:,ifoi,1,:))));
    max_atx(:,2) = squeeze(max(squeeze(atx.perm_k_atx_vox(:,ifoi,2,:))));
    max_dpz(:,1) = squeeze(min(squeeze(dpz.perm_k_dpz_vox(:,ifoi,1,:))));
    max_dpz(:,2) = squeeze(min(squeeze(dpz.perm_k_dpz_vox(:,ifoi,2,:))));
  
    for ivox = 1 : 400
      % PVAL for ATX REST
      p_atx_vox(ivox,ifoi,1) = 1-sum(deg_atx_vox(ivox,ifoi,2)-deg_atx_vox(ivox,ifoi,1)>max_atx(:,1))./nperm;
      % PVAL for ATX TASK
      p_atx_vox(ivox,ifoi,2) = 1-sum(deg_atx_task_vox(ivox,ifoi,2)-deg_atx_task_vox(ivox,ifoi,1)>max_atx(:,2))./nperm;
  %     % PVAL for DPZ REST
      p_dpz_vox(ivox,ifoi,1) = 1-sum(deg_dpz_vox(ivox,ifoi,2)-deg_dpz_vox(ivox,ifoi,1)<max_dpz(:,1))./nperm;
  %     % PVAL for DPZ TASK
      p_dpz_vox(ivox,ifoi,2) = 1-sum(deg_dpz_task_vox(ivox,ifoi,2)-deg_dpz_task_vox(ivox,ifoi,1)<max_dpz(:,2))./nperm;
    end
    
  elseif 1
    % This is already the difference, computed one section earlier  
    max_atx(:,:,1) = squeeze(abs(atx.perm_k_atx_vox(:,ifoi,1,:)));
    max_atx(:,:,2) = squeeze(abs(atx.perm_k_atx_vox(:,ifoi,2,:)));
    max_dpz(:,:,1) = squeeze(abs(dpz.perm_k_dpz_vox(:,ifoi,1,:)));
    max_dpz(:,:,2) = squeeze(abs(dpz.perm_k_dpz_vox(:,ifoi,2,:)));
    
    for ivox = 1 : 400
    % PVAL for ATX REST
    p_atx_vox(ivox,ifoi,1) = 1-sum(abs(deg_atx_vox(ivox,ifoi,2)-deg_atx_vox(ivox,ifoi,1))>max_atx(ivox,:,1))./nperm;
    % PVAL for ATX TASK
    p_atx_vox(ivox,ifoi,2) = 1-sum(abs(deg_atx_task_vox(ivox,ifoi,2)-deg_atx_task_vox(ivox,ifoi,1))>max_atx(ivox,:,2))./nperm;
%     % PVAL for DPZ REST
    p_dpz_vox(ivox,ifoi,1) = 1-sum(abs(deg_dpz_vox(ivox,ifoi,2)-deg_dpz_vox(ivox,ifoi,1))>max_dpz(ivox,:,1))./nperm;
%     % PVAL for DPZ TASK
    p_dpz_vox(ivox,ifoi,2) = 1-sum(abs(deg_dpz_task_vox(ivox,ifoi,2)-deg_dpz_task_vox(ivox,ifoi,1))>max_dpz(ivox,:,2))./nperm;
    end
  end
end

clear perm_k_atx perm_k_dpz perm_k_atx_vox perm_k_dpz_vox

%% PLOT DEGREE ON SURFACE (PLACEBO ONLY)

addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
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

ifoi = 9; icond = 1; %alpha = 0.05;
% fc_avg = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,1,icond,ifoi),1),3),4));
% 
cmap = plasma;

par = nanmean(deg_atx.node_degree(:,ifoi,icond),3);
% par = nanmean(fc_avg(:,ifoi),3);

% cmap      = [cmap(end-15:-1:50,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [5 65];
% para.clim = [0.03 0.08];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn   = sprintf('~/pupmod/plots/pupmod_plot_degree_f%d_c%s_v%d.png',ifoi,num2str(icond),v)
tp_plot_surface(par,para)

%% PLOT DEGREE ON SURFACE (CONTRAST)
% Not 100% clear how to get proper p-values
% Problems: full correction destroys any effects
% being too leniant blows atx effects out of prop

ifoi = 6; alpha = 0.025;

cmap = plasma;

par_atx = deg_atx_task_vox(:,ifoi,2)-deg_atx_task_vox(:,ifoi,1);
par_atx(p_atx_vox(:,ifoi,2)>=alpha) = 0;

cmap      = [cmap(end-15:-1:50,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn   = sprintf('~/pupmod/plots/pupmod_plot_degree_contrast_atx_f%d_v%d.png',ifoi,v)
tp_plot_surface(par_atx,sa_template,para)


ifoi = 7; 

par_dpz = deg_dpz_vox(:,ifoi,2)-deg_dpz_vox(:,ifoi,1);
par_dpz(p_dpz_vox(:,ifoi,1)>=alpha) = 0;

cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(end-15:-1:50,:)];
para      = [];
para.clim = [-0.5 0.5];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn   = sprintf('~/pupmod/plots/pupmod_plot_degree_contrast_dpz_f%d_v%d.png',ifoi,v)
tp_plot_surface(par_dpz,sa_template,para)

