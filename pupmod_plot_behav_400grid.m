%% PLOT EVERYTHING
% pupmod_plot_behav_400grid.m


% Content:

% 1. CORTEX
%---------------
%   1.1. Fraction of altered correlations
%   1.2. Drug effect on FC
%   1.3. Correlations of behavior with FC (per edge)
%   1.4. Task vs. Rest
%---------------
% 2. VTPM (40x40)
%   1.1. Fraction of altered correlations
%   1.2. Drug effect on FC
%   1.3. Correlations of behavior with FC (per edge)
%   1.4. Task vs. Rest
%---------------

%%

clear

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pconn/matlab
addpath ~/pupmod/matlab
v = 12;
fc =pupmod_loadpowcorr(v,0);

para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt = behav;

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);

addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
load sa_meg_template;
grid  = select_chans(sa_meg_template.grid_cortex3000,400); 

% ------------------------------------------
% ------------------------------------------
%    #####    ####     ######   ########## ######## ##      ##
%  ##       ##    ##   ##   ##      ##     ##        ##    ##
% ##       ##      ##  ##   ##      ##     ##         ##  ##
% ##       ##      ##  ######       ##     ######       ##
%  ##       ##    ##   ##   ##      ##     ##        ##   ##
%    #####    ####     ##    ##     ##     ######## ##      ##
%
% ------------------------------------------
% ------------------------------------------
% (1) Correlation of FC with Behavior (during placebo)
% (2) Correlation of FC with Behavior (Drug effect)
% (3) Scatter plot: Mean FC with behavior (significant drug FC only)
% (4) SCATTER PLOT: MEAN FC with behavior (only where behavior correlates with FC)
% (3) Correlations of behavior with FC (per connection)
% (4) Task vs Rest effects (VTPM)
% -----------------------------------------

%% (1) CORRELATION FC WITH BEHAVIOR (PLACEBO): Cortex grid

mask = logical(tril(ones(400,400),-1));
%  = 3;

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ = 1:28; 
icond = 1;

for ipharm = 1 : 3
for ifoi = 1:13
  
  clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
  
  d_fc1      = squeeze(fc(:,:,SUBJ,ipharm,icond,ifoi,1));
  d_fc2      = squeeze(fc(:,:,SUBJ,ipharm,icond,ifoi,2));
  d_behav    = nanmean(behav_cnt(ipharm,SUBJ,:),3);
  d_behav1   = behav_cnt(ipharm,SUBJ,1);
  d_behav2   = behav_cnt(ipharm,SUBJ,2);
  d_fc       = squeeze(nanmean(fc(:,:,SUBJ,ipharm,icond,ifoi,:),7));
  
  % identify and ignore nans
  nan_idx_cnt1   = ~isnan(d_behav1);
  nan_idx_cnt2   = ~isnan(d_behav2);
  nan_idx_meg1   = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
  nan_idx_meg2   = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
  nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
  nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);
  
  d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
  d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
  d_behav     = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);

  [r_cnt,p_cnt]   = tp_corr(d_fc(:,:,nan_idx1),d_behav(:,:,nan_idx1),3);
  [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
  [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);
  
  rr(ifoi) = corr(r_cnt1(mask),r_cnt2(mask))
  figure; set(gcf,'color','w');
  
  subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Counting: Average')
  tp_editplots; set(gca,'FontSize',5)
  subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %d',ifoi))
  tp_editplots; set(gca,'FontSize',5)
  subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %d)',ifoi))
  tp_editplots; set(gca,'FontSize',5);
  subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_fc_corrwithbehav_cortex_cond%d_ipharm%d_f%d_v%d.pdf',icond,ipharm,ifoi,v))
end
end
%% (2) CORRELATION FC WITH BEHAVIOR (DRUG EFFECT): Cortex grid

cmap    = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ    = 1:28; 
% ipharm  = 2;

for ipharm = 2:3
for ifoi = 1:13
  ifoi
  clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
  
  d_fc1      = squeeze(fc(:,:,SUBJ,ipharm,2,ifoi,1))-squeeze(fc(:,:,SUBJ,1,2,ifoi,1));
  d_fc2      = squeeze(fc(:,:,SUBJ,ipharm,2,ifoi,2))-squeeze(fc(:,:,SUBJ,1,2,ifoi,2));
  d_fc       = squeeze(nanmean(fc(:,:,SUBJ,ipharm,2,ifoi,:),7))-squeeze(nanmean(fc(:,:,SUBJ,1,2,ifoi,:),7));

  d_behav    = nanmean(behav_cnt(ipharm,SUBJ,:),3)-nanmean(behav_cnt(1,SUBJ,:),3);
  d_behav1   = behav_cnt(ipharm,SUBJ,1)-behav_cnt(1,SUBJ,1);
  d_behav2   = behav_cnt(ipharm,SUBJ,2)-behav_cnt(1,SUBJ,2);
  
  % identify and ignore nans
  nan_idx_cnt1   = ~isnan(d_behav1);
  nan_idx_cnt2   = ~isnan(d_behav2);
  nan_idx_meg1   = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
  nan_idx_meg2   = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
  nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
  nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);
  
  d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
  d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
  d_behav     = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);

  [r_cnt,p_cnt]   = tp_corr(d_fc(:,:,nan_idx1),d_behav(:,:,nan_idx1),3);
  [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
  [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);
  
  figure; set(gcf,'color','w')
  
  subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Counting: Average')
  tp_editplots; set(gca,'FontSize',5)
  subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %d',ifoi))
  tp_editplots; set(gca,'FontSize',5)
  subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %d)',ifoi))
  tp_editplots; set(gca,'FontSize',5);
  subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_corrwithbehav_cortex_ipharm%d_f%d_v%d.pdf',ipharm,ifoi,v),'-fillpage')
end
end

%% (3) SCATTER PLOT: MEAN FC (only significant connections) with behavior 
foi_range       = unique(round(2.^[1:.5:7]));

close
mask = logical(tril(ones(400,400),-1));
iblock = 1:2;

ipharm = 3;
alpha = 0.05;

figure; set(gcf,'color','w','Position',[100 100 700 800]); 

for ifoi = 1:13

  [h,~,~,s]=ttest(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),7),nanmean(fc(:,:,:,1,2,ifoi,iblock),7),'dim',3,'alpha',alpha);

  idx = h(mask)>0 & s.tstat(mask)>0;

  d = squeeze(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),7)-nanmean(fc(:,:,:,1,2,ifoi,iblock),7));
  for isubj = 1 : 28
    tmp = d(:,:,isubj);
    dd(:,isubj) = tmp(mask);
  end

  m = mean(dd(idx,:));

  d_behav = nanmean(behav_cnt(ipharm,:,iblock),3)-nanmean(behav_cnt(1,:,iblock),3);

  
  subplot(4,4,ifoi);
  scatter(m,d_behav,50,'markerfacecolor','k','markeredgecolor','w');
  tp_editplots; axis square
  lsline; [r,p]=corr(m(:),d_behav(:))
  xlabel('\DeltaFC'); ylabel('\DeltaSwitches');
  title(sprintf('f = %d Hz',foi_range(ifoi)))
  text(double(min(m)),40,sprintf('r=%.3f, p=%.3f',r,p),'fontsize',6)
end

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_corrwithbehav_scatter_pharm%d_v%d.eps',ipharm,v))

%% (4) SCATTER PLOT: MEAN FC with behavior (only where behavior correlates with FC)

ipharm = 1:3;
SUBJ= 1 :28;
icond = 1;
figure; set(gcf,'color','w','Position',[100 100 700 800]); 

for ifoi = 1:13
  
  clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
  
  d_fc1      = squeeze(nanmean(fc(:,:,SUBJ,ipharm,icond,ifoi,1),4));
  d_fc2      = squeeze(nanmean(fc(:,:,SUBJ,ipharm,icond,ifoi,2),4));
  d_fc       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,icond,ifoi,:),7),4));
  
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,:),1),3);
  d_behav1   = nanmean(behav_cnt(ipharm,SUBJ,1),1);
  d_behav2   = nanmean(behav_cnt(ipharm,SUBJ,2),1);
  
  % identify and ignore nans
  nan_idx_cnt1   = ~isnan(d_behav1);
  nan_idx_cnt2   = ~isnan(d_behav2);
  nan_idx_meg1   = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
  nan_idx_meg2   = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
  nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
  nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);
  
  d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
  d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
  d_behav     = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);
  
  [r_cnt,p_cnt]   = tp_corr(d_fc(:,:,nan_idx1),d_behav(:,:,nan_idx1),3);
  [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
  [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);
  
  d_fc       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,2,ifoi,:),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,2,ifoi,:),7),4));
  d_fc1       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,2,ifoi,1),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,2,ifoi,1),7),4));
  d_fc2       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,2,ifoi,2),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,2,ifoi,2),7),4));
  
  for isubj = 1 : 28
   tmp = d_fc(:,:,isubj);
   tmp1 = d_fc1(:,:,isubj);
   tmp2 = d_fc2(:,:,isubj);
   dd(isubj) = mean(tmp((p_cnt<0.05&r_cnt>0)));
   dd1(isubj) = mean(tmp1((p_cnt1<0.05&r_cnt>0)));
   dd2(isubj) = mean(tmp2((p_cnt2<0.05&r_cnt>0)));
  end
  
  dd_b = nanmean(nanmean(behav_cnt(2,SUBJ,:),1),3)-nanmean(nanmean(behav_cnt(1,SUBJ,:),1),3);
  [r,p]=corr(dd(:),dd_b(:));
  
  subplot(4,4,ifoi);
  scatter(dd,dd_b,50,'markerfacecolor','k','markeredgecolor','w');
  tp_editplots; axis square
  lsline;
  xlabel('\DeltaFC'); ylabel('\DeltaSwitches');
  title(sprintf('f = %d Hz',foi_range(ifoi)))
  if ~isnan(r)
    text(double(min(dd)),40,sprintf('r=%.3f, p=%.3f',r,p),'fontsize',6)
  end
end

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_corrwithbehav_scatter_cond%d_pharm%s_v%d.eps',icond,regexprep(num2str(ipharm),' ',''),v))


%% STATISTICS: PLOT SPECTRUM OF CORRELATIONS WITH BEHAVIOR
mask      = logical(tril(ones(400,400),-1));

SUBJ    = 1:28; %SUBJ(13)=[];
iblock  = 1:2;
ipharm  = 3;
cond    = 2;
NPERM   = 10000;
FOI = 6;

% EMPIRICAL
d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav    = permute(repmat(d_behav(:),[1 400 400 13]),[2 3 1 4]);
d_fc       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,:,iblock),7),4));%-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),4));

nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
nan_idx     = nan_idx_meg&nan_idx_beh;

% compute correlation
[r,p] = tp_corr(d_fc(:,:,nan_idx,:),d_behav(:,:,nan_idx,:),3);

for ifoi = FOI
  p_tmp = p(:,:,ifoi); r_tmp = r(:,:,ifoi);
  pos_corr(ifoi) = 100 * (sum(p_tmp(mask)<0.05 & r_tmp(mask)>0) / sum(mask(:)));
  neg_corr(ifoi) = 100 * (sum(p_tmp(mask)<0.05 & r_tmp(mask)<0) / sum(mask(:)));
  pos_corr_vox(:,ifoi) = 100 * (sum(p_tmp<0.05 & r_tmp>0)/size(mask,1));
  neg_corr_vox(:,ifoi) = 100 * (sum(p_tmp<0.05 & r_tmp<0)/size(mask,1));
end

pp_emp = squeeze(p);
rr_emp = squeeze(r);

% PLOT CORRELATION ON CORTICAL SURFACE
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

for ifoi = FOI
  
  para =[];
  para.cmap = cmap;
  % para.cmap = para.cmap(end:-1:1,:);
  para.grid = grid;
  para.dd = 0.75;
  para.clim = [-0.2 0.2];
  para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
  tp_plot_surface(nanmean(rr_emp(:,:,ifoi)),para)
  
% PLOT FRACTION OF SIGNIFICANT CORRELATIONS ON CORTICAL SURFACE
  
%   para =[];
%   para.cmap = plasma;
%   % para.cmap = para.cmap(end:-1:1,:);
%   para.grid = grid;
%   para.dd = 0.75;
%   para.clim = [0 3];
%   para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_masked_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
%   tp_plot_surface(pos_corr_vox(:,ifoi),para)
  
end
%%



%%
mask = logical(tril(ones(400,400),-1));
NPERM = 1000;
clear pp_perm rr_perm r p pos_c*
ipharm = 1:3;
iblock = 1 :2;
cond = 2;

if ~exist(sprintf('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat_pharm%s.mat',regexprep(num2str(ipharm),' ','')))
  for iperm = 1 : NPERM
    
    fn = sprintf('pupmod_plots_behav_400grid_permstat_ipharm%s_iperm%d',regexprep(num2str(ipharm),' ',''),iperm);
    if tp_parallel(fn,'~/pupmod/proc/conn/',1,0)
      continue
    end
    iperm
    order1 = randperm(28);
    order2 = randperm(28);
    
    d_behav    = nanmean(nanmean(behav_cnt(ipharm,order1,iblock),3),1);
    d_behav = permute(repmat(d_behav(:),[1 400 400 13]),[2 3 1 4]);
    d_fc       = squeeze(nanmean(nanmean(fc(:,:,order2,ipharm,cond,:,iblock),7),4));
    
    nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
    nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
    nan_idx     = nan_idx_meg&nan_idx_beh;
    
    % compute correlation
    [r,p]=tp_corr(d_fc(:,:,nan_idx,:),d_behav(:,:,nan_idx,:),3);
    
    for ifoi = 1:13
      p_tmp = p(:,:,:,ifoi);
      r_tmp = r(:,:,:,ifoi);
      out.pos_corr_perm(ifoi) = 100* (sum(p_tmp(mask)<0.05 & r_tmp(mask)>0) / sum(mask(:)));
      out.neg_corr_perm(ifoi) = 100* (sum(p_tmp(mask)<0.05 & r_tmp(mask)<0) / sum(mask(:)));
      out.pos_corr_perm_vox(:,ifoi) = 100*sum(p_tmp<0.05 & r_tmp>0)/size(mask,1);
      out.neg_corr_perm_vox(:,ifoi) = 100*sum(p_tmp<0.05 & r_tmp<0)/size(mask,1);
    end
    
    save(sprintf(['~/pupmod/proc/conn/' fn '.mat']),'out')
    tp_parallel(fn,'~/pupmod/proc/conn/',0,0);
  end
end
cnt = 0;
if ~exist(sprintf('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat_pharm%s.mat',regexprep(num2str(ipharm),' ','')))
  for iperm = 1 : NPERM
    %     iperm
    try
      load(sprintf('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat_pharm%s_iperm%d.mat',regexprep(num2str(ipharm),' ',''),iperm))
      cnt = cnt + 1
      all_pos_corr_perm(cnt,:)=out.pos_corr_perm;
      all_neg_corr_perm(cnt,:)=out.neg_corr_perm;
      all_pos_corr_perm_vox(:,:,cnt)=out.pos_corr_perm_vox;
      all_neg_corr_perm_vox(:,:,cnt)=out.neg_corr_perm_vox;
    catch me
      continue
    end
  end
  save('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat.mat','all_pos_corr_perm','all_neg_corr_perm');
else
  load('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat.mat');
end

error('!')
%%

figure; set(gcf,'color','w');

subplot(2,1,1); hold on

plot(pos_corr,'r');
plot(neg_corr,'b');

tp_editplots;
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
xlabel('Carrier frequency [Hz]'); ylabel('Fraction of sign. correlations [%]')


%%
% diff_rest = nanmean(fc(:,:,:,2,1,ifoi,:),7)-nanmean(fc(:,:,:,1,1,ifoi,:),7);
% diff_task = nanmean(fc(:,:,:,2,2,ifoi,:),7)-nanmean(fc(:,:,:,1,2,ifoi,:),7);
iblock = 1 : 2;
for cond = 2
  for ipharm = 1 : 3
  % ipharm = 3;
  % diff_rest = nanmean(fc(:,:,:,2,1,ifoi,:),7)-nanmean(fc(:,:,:,1,1,ifoi,:),7);
  % diff_task = nanmean(fc(:,:,:,2,2,ifoi,:),7)-nanmean(fc(:,:,:,1,2,ifoi,:),7);
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav    = permute(repmat(d_behav(:),[1 400 13]),[2 1 3]);
  d_fc       = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,:,iblock),2),7),5),4));%-squeeze(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),5),4));

  d_fc = permute(d_fc,[1 3 2]); d_behav = permute(d_behav,[1 3 2]);
  nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
  nan_idx     = nan_idx_meg&nan_idx_beh;

  % compute correlation
  [r,p] = tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);

    for ifoi = 6

      para =[];
      para.cmap = cmap;
      para.grid = grid;
      para.dd = 0.75;
      para.clim = [-0.2 0.2];
      para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_cond%d_b%s_pharm%d_f%d.png',cond,regexprep(num2str(iblock),' ',''),ipharm,ifoi);
      tp_plot_surface(r(:,ifoi),para)

    end
  end
end


%%

for iperm = 1 : 1000
  iperm
  idx = randperm(28);
  
  d_behav    = nanmean(nanmean(behav_cnt(2,idx,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav    = permute(repmat(d_behav(:),[1 400 13]),[2 1 3]);
  d_fc       = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,idx,2,cond,:,iblock),2),7),5),4));%-squeeze(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),5),4));

  d_fc = permute(d_fc,[1 3 2]); d_behav = permute(d_behav,[1 3 2]);
  nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
  nan_idx     = nan_idx_meg&nan_idx_beh;

  [r2,p] = tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);

  d_behav    = nanmean(nanmean(behav_cnt(1,idx,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav    = permute(repmat(d_behav(:),[1 400 13]),[2 1 3]);
  d_fc       = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,idx,1,cond,:,iblock),2),7),5),4));%-squeeze(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),5),4));

  d_fc = permute(d_fc,[1 3 2]); d_behav = permute(d_behav,[1 3 2]);
  nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
  nan_idx     = nan_idx_meg&nan_idx_beh;

  [r1,p] = tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);

  d(:,:,iperm) = r2-r1;
end


% EMPRICIAL RESULT
d_behav    = nanmean(nanmean(behav_cnt(2,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav    = permute(repmat(d_behav(:),[1 400 13]),[2 1 3]);
d_fc       = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,SUBJ,2,cond,:,iblock),2),7),5),4));%-squeeze(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),5),4));

d_fc = permute(d_fc,[1 3 2]); d_behav = permute(d_behav,[1 3 2]);
nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
nan_idx     = nan_idx_meg&nan_idx_beh;

[r2,p] = tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);

d_behav    = nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav    = permute(repmat(d_behav(:),[1 400 13]),[2 1 3]);
d_fc       = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),2),7),5),4));%-squeeze(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,:,iblock),7),5),4));

d_fc = permute(d_fc,[1 3 2]); d_behav = permute(d_behav,[1 3 2]);
nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
nan_idx     = nan_idx_meg&nan_idx_beh;

[r1,p] = tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);


% PLOT
p =1-(sum((r2(:,6)-r1(:,6))>squeeze(d(:,6,:)),2)/size(d,3));

para =[];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.35 0.35];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_cond%d_b%s_pharm%d_f%d.png',cond,regexprep(num2str(iblock),' ',''),ipharm,ifoi);
tp_plot_surface((r2(:,6)-r1(:,6)).*(p<0.1),para)

%%
% ------------------------------------------
% ------------------------------------------
%
% ##       ## ######### ######   ##       ##
%  ##     ##     ##     ##   ##  ###     ###
%   ##   ##      ##     ##   ##  ## #   # ##
%    ## ##       ##     #####    ##  # #  ##
%     ##         ##     ##       ##   #   ##
%
% ------------------------------------------
% ------------------------------------------
% (1) Fraction of altered correlations (VTPM full)
% (2) Drug effects on FC (VTPM full)
% (3) Correlations of behavior with FC (per connection)
% (4) Task vs Rest effects (VTPM)
% ------------------------------------------

% load data
fc_vtpm = pupmod_loadpowcorr(20,0);

tmp = tp_create_grid('vtpm');
idx = 1:46;
idx = ~ismember(idx,[21 22 23 44 45 46]);
reg = tmp.tissuelabel_4mm(idx);

%% (1) FRACTION OF ALTERED CORRELATION IN VTPM
% -------------------------------

mask = logical(tril(ones(40,40),-1));

for ipharm = 2:3
  for cond = 1 : 2
    for ifoi = 1 : 13
      
      [h,~,~,s]= ttest(nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,cond,ifoi,:),6),7),nanmean(nanmean(fc_vtpm(idx,idx,:,1,cond,ifoi,:),6),7),'dim',3);
      
      pos(ifoi,cond,ipharm-1) = 100*sum((h(mask)>0 & s.tstat(mask)>0))/sum(mask(:));
      neg(ifoi,cond,ipharm-1) = 100*sum((h(mask)>0 & s.tstat(mask)<0))/sum(mask(:));
      
    end
  end
end

figure; set(gcf,'color','w');

subplot(3,2,1); hold on
plot(pos(:,1,1),'r');
plot(neg(:,1,1),'b');
axis([0 14 0 100])

set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
tp_editplots; title('Rest | Atx')
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of\naltered corr.'))

subplot(3,2,3); hold on
plot(pos(:,2,1),'r');
plot(neg(:,2,1),'b');
axis([0 14 0 100])

set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
tp_editplots; title('Rest | Dpz')
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of\naltered corr.'))

subplot(3,2,2); hold on
plot(pos(:,1,2),'r');
plot(neg(:,1,2),'b');
axis([0 14 0 100])

set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
tp_editplots; title('Task | Atx')
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of\naltered corr.'))

subplot(3,2,4); hold on
plot(pos(:,2,2),'r');
plot(neg(:,2,2),'b');
axis([0 14 0 100])

set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
tp_editplots;  title('Task | Dpz')
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of\naltered corr.'))

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_plot_behav_alteredcorr_vtpm.eps'))

%% (2) PLOT DRUG EFFECTS ON FC IN VTPM
% ----------------------------------------
% During Task (icond = 2) or Rest (icond = 1)
% for ipharm = 2 (Atx vs Pbo) and ipharm = 3 (Dpz vs Pbo)

mask = logical(tril(ones(40,40),-1)); 
idx = 1:46; idx = ~ismember(idx,[21 22 23 44 45 46]);

for ipharm = 2 : 3
  for ifoi = 1 : 13
    figure; set(gcf,'color','w');
    
    for icond = 1 : 2
      
      [h,p]= ttest(nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,:),6),7),nanmean(nanmean(fc_vtpm(idx,idx,:,1,icond,ifoi,:),6),7),'dim',3);
      h = p<(fdr1(p(mask),0.05));
      h = p<0.05;
      
      d = nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,:),3),7)-nanmean(nanmean(fc_vtpm(idx,idx,:,1,icond,ifoi,:),3),7);
      
      d = nanmean(d,3);
      clim = [-max([abs(min(d(:))), abs(max(d(:)))]) max([abs(min(d(:))), abs(max(d(:)))])];
      
      subplot(2,2,1+(icond-1));
      imagesc(d,clim)
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square; title(sprintf('Drug%d vs Drug1: Cond%d',ipharm,cond))
      tp_colorbar()
      
      subplot(2,2,3+(icond-1));
      imagesc(d.*h,clim)
      
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square
      title(sprintf('Masked at P=0.05 (Cond%d)',cond))
      
      colormap(cmap)
      
      print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_plot_drugonfc_vtpm_full_ipharm%d_f%d.eps',ipharm,ifoi))
      
    end
  end
end

% ----------------------------
%% (3) PLOT CORRELATIONS BEHAVIOR (PBO) VS FC (PBO)
% ----------------------------
% Correlate switch counts and FC across subjects, for each connections
% Plots unthresholded and thresholded correlation matrices (P=0.05)

mask = logical(tril(ones(40,40),-1)); 
idx = 1:46; idx = ~ismember(idx,[21 22 23 44 45 46]);
foi_range       = unique(round(2.^[1:.5:7]));
clim = [-0.5 0.5];

icond = 2;

for ipharm = 1 : 3
  
  for ifoi = 1 : 13
    figure; set(gcf,'color','w');
    for iblock = 1 : 2
        
      d_behav =  nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1);
      d = nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,iblock),7),4);
      
      nanidx = ~isnan(d_behav') & ~isnan(squeeze(d(1,2,:)));
      
      clear r p
      for i = 1 : 40
        for j = 1 : 40
          [r(i,j) p(i,j)] = corr(squeeze(d(i,j,nanidx)),d_behav(nanidx)');
        end
      end
      
      subplot(3,2,2*(iblock-1)+1);
      imagesc(r,clim)
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
      tp_colorbar()
      
      subplot(3,2,2*(iblock-1)+2);
      imagesc(r.*(p<0.05),clim)
      
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square
      title(sprintf('Masked at P=0.05 (Block%d)',iblock))
      
      colormap(cmap)
      
    end
    
    d_behav =  nanmean(nanmean(behav_cnt(ipharm,:,:),3),1);
    d = nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,:),7),4);
    
    clear r p
    for i = 1 : 40
      for j = 1 : 40
        [r(i,j) p(i,j)] = corr(squeeze(d(i,j,:)),d_behav(:));
      end
    end
    
    subplot(3,2,5);
    imagesc(r,clim)
    tp_editplots;
    set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
    set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
    set(gca,'FontSize',5);
    axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
    tp_colorbar()
    
    subplot(3,2,6);
    imagesc(r.*(p<0.05),clim)
    tp_editplots;
    set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
    set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
    set(gca,'FontSize',5);
    axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
    tp_colorbar()
    
    %   drawnow
    print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_fc_corrwithbehav_vtpm_full_pharm_%d_f%d.eps',ipharm,ifoi))
    
  end
end


% --------------------------------
%% (4) CORRELATE DRUG EFFECTS AND BEHAVIORAL EFFECTS (VTPM)

mask = logical(tril(ones(40,40),-1)); 
idx = 1:46; idx = ~ismember(idx,[21 22 23 44 45 46]);
foi_range       = unique(round(2.^[1:.5:7]));
clim = [-0.5 0.5];

icond = 2;

for ipharm = 2 : 3
  
  for ifoi = 1 : 13
    figure; set(gcf,'color','w');
    for iblock = 1 : 2
        
      d_behav =  nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1)-nanmean(nanmean(behav_cnt(1,:,iblock),3),1);
      d = nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,iblock),7),4)-nanmean(nanmean(fc_vtpm(idx,idx,:,1,icond,ifoi,iblock),7),4);
      
      nanidx = ~isnan(d_behav') & ~isnan(squeeze(d(1,2,:)));
      
      clear r p
      for i = 1 : 40
        for j = 1 : 40
          [r(i,j) p(i,j)] = corr(squeeze(d(i,j,nanidx)),d_behav(nanidx)');
        end
      end
      
      subplot(3,2,2*(iblock-1)+1);
      imagesc(r,clim)
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
      tp_colorbar()
      
      subplot(3,2,2*(iblock-1)+2);
      imagesc(r.*(p<0.05),clim)
      
      tp_editplots;
      set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
      set(gca,'FontSize',5);
      axis square
      title(sprintf('Masked at P=0.05 (Block%d)',iblock))
      
      colormap(cmap)
      
    end
    
    d_behav =  nanmean(nanmean(behav_cnt(ipharm,:,:),3),1);
    d = nanmean(nanmean(fc_vtpm(idx,idx,:,ipharm,icond,ifoi,:),7),4);
    
    clear r p
    for i = 1 : 40
      for j = 1 : 40
        [r(i,j) p(i,j)] = corr(squeeze(d(i,j,:)),d_behav(:));
      end
    end
    
    subplot(3,2,5);
    imagesc(r,clim)
    tp_editplots;
    set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
    set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
    set(gca,'FontSize',5);
    axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
    tp_colorbar()
    
    subplot(3,2,6);
    imagesc(r.*(p<0.05),clim)
    tp_editplots;
    set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
    set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
    set(gca,'FontSize',5);
    axis square; title(sprintf('Corr FC-Behav (Pbo): f = %d Hz',foi_range(ifoi)))
    tp_colorbar()
    
    %   drawnow
    print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_corrwithbehav_vtpm_full_pharm_%d_f%d.eps',ipharm,ifoi))
    
  end
end



%% COARSER REGIONS
% labels = {'V1';'V2-4';'V3ab';'IPS01';'IPS23';'LO';'VO';'PHC';'MT'};
% region{1} = [1 2]; region{2} = [3 4 5 6 7]; region{3} = [15 16]; region{4} = [17 18];
% region{5} = [19 20]; region{6} = [13 14]; region{7} = [8 9]; region{8} = [10 11]; region{9} = [12];
%
% fc_vtpm_lr = squeeze((fc_vtpm(1:20,1:20,:,:,:,:,:)+fc_vtpm(24:43,24:43,:,:,:,:,:))./2);
% clear tmp fc
% %
% % % collapse across regions defined in 'labels'
% % for i = 1 : size(region,2)
% %   for j = 1 : size(region,2)
% %     if i == j; fc(i,j,:,:,:) = nan(28,3,2); continue; end
% %     clear tmp
% %     for ii = 1 : length(region{i})
% %       tmp(ii,:,:,:) = squeeze(mean(cleandat_lr(region{i}(ii),region{j},:,:,:),2));
% %     end
% %     fc(i,j,:,:,:) = squeeze(mean(tmp,1));
% %   end
% % end
%
% cond = 2;
% ipharm = 2;
% iblock = 1:2;
%
% % diff_task = nanmean(fc(:,:,:,2,2,ifoi,:),7)-nanmean(fc(:,:,:,1,2,ifoi,:),7);
% d_behav    = nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1)-nanmean(nanmean(behav_cnt(1,:,iblock),3),1);
% d_behav    = permute(repmat(d_behav(:),[1 20 20 13]),[2 3 1 4]);
% d_fc       = squeeze(nanmean(nanmean(nanmean(fc_vtpm_lr(:,:,:,ipharm,cond,:,iblock),7),5),4))-squeeze(nanmean(nanmean(nanmean(fc_vtpm_lr(:,:,:,1,cond,:,iblock),7),5),4));
%
% nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:,1)));
% nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:,1)));
% nan_idx     = nan_idx_meg&nan_idx_beh;
%
% % compute correlation
% [r,p] = tp_corr(d_fc(:,:,nan_idx,:),d_behav(:,:,nan_idx,:),3);
% %%
% clim = [-0.25 0.25];
% ifoi = 2;
% figure; set(gcf,'color','w');
% subplot(1,2,1);
% imagesc(r(:,:,:,ifoi),clim)
% tp_editplots;
% set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
% set(gca,'FontSize',7);
%
% subplot(1,2,2);
% imagesc(r(:,:,:,ifoi).*(p(:,:,:,ifoi)<0.05),clim)
%
% tp_editplots;
% set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
% set(gca,'FontSize',7);
%
% colormap(cmap)


% ------------------------------------------
% ------------------------------------------
%
% ##       ## ######### ######   ##       ##        ##       ######
%  ##     ##     ##     ##   ##  ###     ###        ##       ##   ##
%   ##   ##      ##     ##   ##  ## #   # ##  ####  ##       ######
%    ## ##       ##     #####    ##  # #  ##        ##       ##  ##
%     ##         ##     ##       ##   #   ##        #######  ##   ##
%
% ------------------------------------------
% ------------------------------------------

