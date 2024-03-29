edit  %% PLOT EVERYTHING
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

v  = 23;
cleandat = pupmod_loadpowcorr(v,SUBJLIST,0);
% fc = cat(4,fc(:,:,:,:,:,:,1),fc(:,:,:,:,:,:,2));
% fc = (fc-nanmean(fc,4))./nanstd(fc,[],4);
% fc = cat(7,fc(:,:,:,1:3,:,:),fc(:,:,:,4:6,:,:));
%
para.str_behav  = 'count';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt       = behav;

tmp = reshape(behav_cnt,[3*2 28]);
tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
behav_cnt = reshape(tmp,[3 28 2]);

para.str_behav  = 'numb_switches';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn      = behav;
behav_bttn      = permute(behav_bttn,[2 1 3]);
%
tmp = reshape(behav_bttn,[3*2 28]);
tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
behav_bttn = reshape(tmp,[3 28 2]);

behav_cnt(isnan(behav_cnt))   = behav_bttn(isnan(behav_cnt));
behav_bttn(isnan(behav_bttn)) = behav_cnt(isnan(behav_bttn));

behav_cnt = (behav_cnt+behav_bttn)./2;

addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
load sa_meg_template.mat;
grid  = select_chans(sa_meg_template.grid_cortex3000,400);
clear sa_meg_template

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
clear rr
mask = logical(tril(ones(400,400),-1));

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ = 1:28;
icond = 2;
foi_range = 2.^(1:0.25:7);
pos = [0.63 0.63 19.73 28.43];
ipharm = 1;

for ifoi = 1:21
  
  clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
  
%   d_fc1      = squeeze(nanmean(cleandat(:,:,SUBJ,ipharm,icond,ifoi,1),4));
%   d_fc2      = squeeze(nanmean(cleandat(:,:,SUBJ,ipharm,icond,ifoi,2),4));
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,:),3),1);
%   d_behav1   = nanmean(behav_cnt(ipharm,SUBJ,1),1);
%   d_behav2   = nanmean(behav_cnt(ipharm,SUBJ,2),1);
  d_fc       = squeeze(nanmean(cleandat(:,:,SUBJ,ipharm,icond,ifoi,:),4));
  
  % identify and ignore nans
%   nan_idx_cnt1   = ~isnan(d_behav1);
%   nan_idx_cnt2   = ~isnan(d_behav2);
%   nan_idx_meg1   = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
%   nan_idx_meg2   = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
%   nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
%   nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);
  
%   d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
%   d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
  d_behav     = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);
  
  [r_cnt,p_cnt]   = tp_corr(d_fc,d_behav,3);
%   [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
%   [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);

%   rr(ifoi) = corr(r_cnt1(mask),r_cnt2(mask))
  figure; set(gcf,'color','w');
  subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title(sprintf('Counting: Average (f=%d',ifoi))
  tp_editplots; set(gca,'FontSize',5)
  subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
  colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
  
%   subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %.2f Hz)',foi_range(ifoi)))
%   tp_editplots; set(gca,'FontSize',5)
%   subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
%   colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
%   
%   subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %.2f Hz)',foi_range(ifoi)))
%   tp_editplots; set(gca,'FontSize',5);
%   subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
%   colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
%   
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_fc_corrwithbehav_cortex_cond%d_ipharm%s_f%d_v%d.pdf',icond,regexprep(num2str(ipharm),' ',''),ifoi,v),'-fillpage')
  close
  
end
%% (2) CORRELATION FC WITH BEHAVIOR (DRUG EFFECT): Cortex grid

cmap    = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ    = 1:28;

cond = 2;
for ipharm = 2
  for ifoi = 1:25
    ifoi
    clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
    
%     d_fc1      = squeeze(fc(:,:,SUBJ,ipharm,cond,ifoi,1))-squeeze(fc(:,:,SUBJ,1,cond,ifoi,1));
%     d_fc2      = squeeze(fc(:,:,SUBJ,ipharm,cond,ifoi,2))-squeeze(fc(:,:,SUBJ,1,cond,ifoi,2));
    d_fc       = squeeze(cleandat(:,:,SUBJ,ipharm,cond,ifoi))-squeeze(cleandat(:,:,SUBJ,1,cond,ifoi));
    
    d_behav    = nanmean(behav_cnt(ipharm,SUBJ,:),3)-nanmean(behav_cnt(1,SUBJ,:),3);
%     d_behav1   = behav_cnt(ipharm,SUBJ,1)-behav_cnt(1,SUBJ,1);
%     d_behav2   = behav_cnt(ipharm,SUBJ,2)-behav_cnt(1,SUBJ,2);
    
    % identify and ignore nans
%     nan_idx_cnt1   = ~isnan(d_behav1);
%     nan_idx_cnt2   = ~isnan(d_behav2);
%     nan_idx_meg1   = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
%     nan_idx_meg2   = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
%     nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
%     nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);
    
%     d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
%     d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
    d_behav     = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);
    
    [r_cnt,p_cnt]   = tp_corr(d_fc,d_behav,3);
%     [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
%     [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);
    
    figure; set(gcf,'color','w')
    
    subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Counting: Average')
    tp_editplots; set(gca,'FontSize',5)
    subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
    colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
    
%     subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %.2f)',foi_range(ifoi)))
%     tp_editplots; set(gca,'FontSize',5)
%     subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
%     colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
%     
%     subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %.2f)',foi_range(ifoi)))
%     tp_editplots; set(gca,'FontSize',5);
%     subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
%     colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');
    
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_corrwithbehav_cortex_cond%d_ipharm%d_f%d_v%d.pdf',cond,ipharm,ifoi,v),'-fillpage')
    close
  end
end

%% (3) SCATTER PLOT: MEAN FC (only significant connections) with behavior
foi_range       = 2.^[1:.25:7];

close
mask    = logical(tril(ones(400,400),-1));
iblock  = 1:2;
% ipharm = 1;
cond    = 2;
alpha   = 0.05;

for ipharm = 2 :3
  figure; set(gcf,'color','w','Position',[100 100 700 800]);
  
  for ifoi = 1:25
    
    [h,~,~,s]=ttest(nanmean(cleandat(:,:,:,ipharm,cond,ifoi,iblock),7),nanmean(cleandat(:,:,:,1,cond,ifoi,iblock),7),'dim',3,'alpha',alpha);
    
    idx = h(mask)>0 & s.tstat(mask)>0;
    
    d = squeeze(nanmean(cleandat(:,:,:,ipharm,cond,ifoi,iblock),7)-nanmean(cleandat(:,:,:,1,cond,ifoi,iblock),7));
    for isubj = 1 : 28
      tmp = d(:,:,isubj);
      dd(:,isubj) = tmp(mask);
    end
    
    d_fc1      = squeeze(nanmean(cleandat(:,:,SUBJ,ipharm,icond,ifoi,1),4));
    d_fc2      = squeeze(nanmean(cleandat(:,:,SUBJ,ipharm,icond,ifoi,2),4));
    m = mean(dd(idx,:));
    
    d_behav = nanmean(behav_cnt(ipharm,:,iblock),3)-nanmean(behav_cnt(1,:,iblock),3);
    
    subplot(5,5,ifoi);
    scatter(m,d_behav,50,'markerfacecolor','k','markeredgecolor','w');
    tp_editplots; axis square
    lsline; [r,p]=corr(m(:),d_behav(:))
    xlabel('\DeltaFC'); ylabel('\DeltaSwitches');
    title(sprintf('f = %.2f Hz',foi_range(ifoi)))
    text(double(min(m)),-2,sprintf('r=%.3f, p=%.3f',r,p),'fontsize',6)
    axis([double(min(m)) double(max(m)) -2.5 2.5])
  end
  
  print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_signconn_corrwithbehav_scatter_cond%d_pharm%d_v%d.eps',cond,ipharm,v))
end
%% (4) SCATTER PLOT: MEAN FC with behavior (only where behavior correlates with FC)

clear dd dd1 dd2
ipharm  = 1:3;
SUBJ    = 1:28;
cond    = 2;

figure; set(gcf,'color','w','Position',[100 100 700 800]);

for ifoi = 1:25
  
  clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2
  
  d_fc1      = squeeze(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,1),4));
  d_fc2      = squeeze(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,2),4));
  d_fc       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,:),7),4));
  
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,:),1),3);
  d_behav1   = nanmean(behav_cnt(ipharm,SUBJ,1),1);
  d_behav2   = nanmean(behav_cnt(ipharm,SUBJ,2),1);
  
  % identify and ignore nans
  nan_idx_cnt1	= ~isnan(d_behav1);
  nan_idx_cnt2  = ~isnan(d_behav2);
  nan_idx_meg1  = ~any(isnan(squeeze(d_fc1(1,2,:))),2);
  nan_idx_meg2  = ~any(isnan(squeeze(d_fc2(1,2,:))),2);
  nan_idx1      = nan_idx_cnt1(:)&nan_idx_meg1(:);
  nan_idx2      = nan_idx_cnt2(:)&nan_idx_meg2(:);
  
  d_behav1	= permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
  d_behav2	= permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
  d_behav 	= permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);
  
  [r_cnt,p_cnt]   = tp_corr(d_fc(:,:,nan_idx1),d_behav(:,:,nan_idx1),3);
  [r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
  [r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx2),d_behav2(:,:,nan_idx2),3);
  
  d_fc 	= squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,:),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,cond,ifoi,:),7),4));
  d_fc1	= squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,1),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,cond,ifoi,1),7),4));
  d_fc2	= squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,cond,ifoi,2),7),4))-squeeze(nanmean(nanmean(fc(:,:,SUBJ,1,cond,ifoi,2),7),4));
  
  for isubj = 1 : 28
    tmp        = d_fc(:,:,isubj);
    tmp1       = d_fc1(:,:,isubj);
    tmp2       = d_fc2(:,:,isubj);
    dd(isubj)  = mean(tmp((p_cnt<0.05&r_cnt>0)));
    dd1(isubj) = mean(tmp1((p_cnt1<0.05&r_cnt>0)));
    dd2(isubj) = mean(tmp2((p_cnt2<0.05&r_cnt>0)));
  end
  
  dd_b = nanmean(nanmean(behav_cnt(ipharm,SUBJ,:),1),3)-nanmean(nanmean(behav_cnt(1,SUBJ,:),1),3);
  [r,p]=corr(dd(:),dd_b(:))
  
  subplot(5,5,ifoi);
  scatter(dd,dd_b,50,'markerfacecolor','k','markeredgecolor','w');
  tp_editplots; axis square
  lsline;
  xlabel('\DeltaFC'); ylabel('\DeltaSwitches');
  title(sprintf('f = %.2f Hz',foi_range(ifoi)))
  if ~isnan(r)
    text(double(min(dd)),20,sprintf('r=%.3f, p=%.3f',r,p),'fontsize',6)
  end
  axis([double(min(m)) double(max(m)) -2.5 2.5])
end

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_drugeffectonfc_signcorrbehav_corrwithbehav_scatter_cond%d_pharm%s_v%d.eps',icond,regexprep(num2str(ipharm),' ',''),v))

%% SURFACE PLOTS: VOXEL-WISE CORRELATIONS BETWEEN FC AND BEHAVIOR
% ----------------
mask    = logical(tril(ones(400,400),-1));
SUBJ    = 1:28; %SUBJ(13)=[];
iblock  = 1:2;
ipharm  = 1;
cond    = 2;
% FOI     = 25;

for FOI = 12
  
  % EMPIRICAL
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav    = permute(repmat(d_behav,[400 1]),[1 2]);
  d_fc       = squeeze(nanmean(nanmean(nanmean(cleandat(:,:,SUBJ,ipharm,cond,FOI),6),4)));%-squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,cond,FOI,iblock),6),7),4)));
  d_fc       = permute(d_fc,[1 3 2]);
  
  nan_idx_meg = ~isnan(squeeze(d_fc(1,:)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,:)));
  nan_idx     = nan_idx_meg&nan_idx_beh;
  
  % compute correlation
  [r,p] = tp_corr(d_fc(:,nan_idx),d_behav(:,nan_idx),2);
  
  % PLOT CORRELATION ON CORTICAL SURFACE
  cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
  
  para =[];
  para.cmap = cmap;
  para.grid = grid;
  para.dd = 0.75;
  para.clim = [-0.4 0.4];
  para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_within_f%s_b%s_pharm%s_v%d.png',regexprep(num2str(ifoi),' ',''),regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
  tp_plot_surface(r,para)
  
  % Masked
  para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_within_masked_f%s_b%s_pharm%s_v%d.png',regexprep(num2str(ifoi),' ',''),regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
  tp_plot_surface(r.*(p<0.05),para)
  
  close all
end
%% SURFACE PLOTS: VOXEL-WISE CORRELATIONS BETWEEN FC AND BEHAVIOR (DIFFERENCES)
% ----------------
mask    = logical(tril(ones(400,400),-1));
SUBJ    = 1:28; %SUBJ(13)=[];
iblock  = 1:2;
ipharm  = 2;
cond    = 2;
FOI     = 24;

% for FOI = FOIS
  % EMPIRICAL
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1)-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav    = permute(repmat(d_behav,[400 1]),[1 2]);
  d_fc       = squeeze(nanmean(nanmean(nanmean(cleandat(:,:,SUBJ,ipharm,cond,FOI),6),4)))-squeeze(nanmean(nanmean(nanmean(cleandat(:,:,SUBJ,1,cond,FOI),6),4)));
  d_fc       = permute(d_fc,[1 3 2]);
  
  nan_idx_meg = ~isnan(squeeze(d_fc(1,:)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,:)));
  nan_idx     = nan_idx_meg&nan_idx_beh;
  
  % compute correlation
  [r,p] = tp_corr(d_fc(:,nan_idx),d_behav(:,nan_idx),2);
  
  % PLOT CORRELATION ON CORTICAL SURFACE
  cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
  
  para =[];
  para.cmap = cmap;
  para.grid = grid;
  para.dd = 0.75;
  para.clim = [-0.4 0.4];
  para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_differences_f%s_b%s_pharm%s_v%d.png',regexprep(num2str(FOI),' ',''),regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
  tp_plot_surface(r,para)
  
  % Masked
  para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_differences_masked_f%s_b%s_pharm%s_v%d.png',regexprep(num2str(FOI),' ',''),regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
  tp_plot_surface(r.*(p<0.05),para)
  
%   close all
% end

%% SPECTRUM OF CORRELATIONS BETWEEN FC AND BEHAVIOR
% ---------------------------------

% load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
% lab=tp_aal2yeo(grid);

% idx = lab==1;
idx=ones(400,1)==1;
mask    = logical(triu(ones(sum(idx),sum(idx)),1));

SUBJ    = 1:28;
iblock  = 1:2;
ipharm  = 2;
cond    = 2;
NPERM   = 10000;
FOI     = 1:25;

d_behav    = nanmean(nanmean(behav_bttn(ipharm,SUBJ,iblock),3),1)-nanmean(nanmean(behav_bttn(1,SUBJ,iblock),3),1);
d_behav2   = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1)-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav = (d_behav+d_behav2)./2;
d_behav(27) = d_behav2(27);

SUBJ = d_behav2>0;

% EMPIRICAL
d_behav    = nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav    = permute(repmat(d_behav(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
d_fc       = squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
d_fc       = permute(d_fc,[1 2 4 3]);

nan_idx_meg = ~isnan(squeeze(d_fc(1,2,1,:)));
nan_idx_beh = ~isnan(squeeze(d_behav(1,2,1,:)));
nan_idx     = nan_idx_meg&nan_idx_beh;

% compute correlation
[r,p] = tp_corr(d_fc(:,:,:,nan_idx),d_behav(:,:,:,nan_idx),4);
clear pos_corr_vox neg_corr_vox

for ifoi = FOI
  p_tmp = p(:,:,ifoi); r_tmp = r(:,:,ifoi);
  pos_corr(ifoi) = 100 * (sum(p_tmp(mask)<0.05 & r_tmp(mask)>0) / sum(mask(:)));
  neg_corr(ifoi) = 100 * (sum(p_tmp(mask)<0.05 & r_tmp(mask)<0) / sum(mask(:)));
  pos_corr_vox(:,ifoi) = 100 * (sum(p_tmp<0.05 & r_tmp>0)/size(mask,1));
  neg_corr_vox(:,ifoi) = 100 * (sum(p_tmp<0.05 & r_tmp<0)/size(mask,1));
  
  
  tmp = r(:,:,ifoi);
  tmp = tmp((r(:,:,ifoi)>0));
  r_pos(ifoi) = mean(tmp);
  
  tmp = r(:,:,ifoi);
  tmp = tmp((r(:,:,ifoi)<0));
  r_neg(ifoi) = mean(tmp);
end

% for isubj = 1 : 28
%   for ifoi = 1 : 25
%     tmp = d_fc(:,:,ifoi,isubj);
%     all(isubj,ifoi) = mean(tmp(p(:,:,ifoi)<0.05));
%   end
% end


figure; set(gcf,'color','w')
subplot(3,2,1); hold on
plot(pos_corr,'r');
plot(neg_corr,'b');
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Fraction of sign. correlations [%]'); tp_editplots
tp_editplots
axis([1 25 0 30])

subplot(3,2,2); hold on
plot(r_pos,'r');
plot(abs(r_neg),'b');
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Fraction of sign. correlations [%]'); tp_editplots
tp_editplots
axis([1 25 0 0.6])
% for ifoi = FOI

%%
iblock = 1:2;
% idx = lab==1
idx=ones(400,1)==1;
mask    = logical(triu(ones(sum(idx),sum(idx)),1));

% d_behav    = nanmean(nanmean(behav_bttn(2,:,iblock),3),1)-nanmean(nanmean(behav_bttn(1,:,iblock),3),1);
% d_behav2   = nanmean(nanmean(behav_cnt(2,:,iblock),3),1)-nanmean(nanmean(behav_cnt(1,:,iblock),3),1);
d_behav2   = nanmean(nanmean(behav_bttn(2,:,iblock),3),1)-nanmean(nanmean(behav_bttn(1,:,iblock),3),1);

% d_behav = (d_behav+d_behav2)./2;
% d_behav(27) = d_behav2(27);
% d_behav2 = d_behav;

[i,k]=sort(d_behav2,'descend')

cond = 2;

SUBJ1 = k(1:14);
SUBJ2 = k(15:28);
% SUBJ3 = k(20:28);

% SUBJ1 = d_behav2>0;
% SUBJ2 = d_behav2<0;

ipharm = 1;

d_behav1    = nanmean(nanmean(behav_cnt(ipharm,SUBJ1,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav1    = permute(repmat(d_behav1(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
d_fc1       = squeeze(nanmean(cleandat(idx,idx,SUBJ1,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
d_fc1       = permute(d_fc1,[1 2 4 3]);

d_behav2    = nanmean(nanmean(behav_cnt(ipharm,SUBJ2,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
d_behav2    = permute(repmat(d_behav2(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
d_fc2       = squeeze(nanmean(cleandat(idx,idx,SUBJ2,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
d_fc2       = permute(d_fc2,[1 2 4 3]);

% d_behav3    = nanmean(nanmean(behav_cnt(ipharm,SUBJ3,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
% d_behav3    = permute(repmat(d_behav3(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
% d_fc3       = squ% d_behav3    = nanmean(nanmean(behav_cnt(ipharm,SUBJ3,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
% d_behav3    = permute(repmat(d_behav3(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
% d_fc3       = squeeze(nanmean(cleandat(idx,idx,SUBJ3,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
% d_fc3       = permute(d_fc3,[1 2 4 3]);
% eeze(nanmean(cleandat(idx,idx,SUBJ3,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
% d_fc3       = permute(d_fc3,[1 2 4 3]);

% d_behav4    = nanmean(nanmean(behav_cnt(ipharm,SUBJ4,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
% d_behav4    = permute(repmat(d_behav4(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
% d_fc4       = squeeze(nanmean(cleandat(idx,idx,SUBJ4,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
% d_fc4       = permute(d_fc4,[1 2 4 3]);



[r1,p1] = tp_corr(d_fc1,d_behav1,4);
[r_all1,p_all1] = tp_corr(squeeze(nanmean(nanmean(d_fc1,1),2)),squeeze(nanmean(nanmean(d_behav1,1),2)),2);

[r2,p2] = tp_corr(d_fc2,d_behav2,4);
[r_all2,p_all2] = tp_corr(squeeze(nanmean(nanmean(d_fc2,1),2)),squeeze(nanmean(nanmean(d_behav2,1),2)),2);

% [r3,p3] = tp_corr(d_fc3,d_behav3,4);
% [r_all3,p_all3] = tp_corr(squeeze(nanmean(nanmean(d_fc3,1),2)),squeeze(nanmean(nanmean(d_behav3,1),2)),2);

% [r4,p4] = tp_corr(d_fc4,d_behav4,4);
% [r_all4,p_all4] = tp_corr(squeeze(nanmean(nanmean(d_fc4,1),2)),squeeze(nanmean(nanmean(d_behav4,1),2)),2);

% 
for ifoi = 1:25
  tmp = r1(:,:,ifoi);
  tmp = tmp((r1(:,:,ifoi)>0));
  r_pos1(ifoi) = tanh(mean(atanh(tmp)));
  k(20:28)
  tmp = r1(:,:,ifoi);
  tmp = tmp((r1(:,:,ifoi)<0));
  r_neg1(ifoi) = tanh(mean(atanh(tmp)));
  
  tmp = r2(:,:,ifoi);
  tmp = tmp((r2(:,:,ifoi)>0));
  r_pos2(ifoi) = tanh(mean(atanh(tmp)));
  
  tmp = r2(:,:,ifoi);
  tmp = tmp((r2(:,:,ifoi)<0));
  r_neg2(ifoi) = tanh(mean(atanh(tmp)));
  
 
end


figure; set(gcf,'color','w')
subplot(3,2,1); hold on; title('Edge-wise correlation')
plot(squeeze(nanmean(nanmean(r1,1),2)),'color',[0.5 0.2 0])
plot(squeeze(nanmean(nanmean(r2,1),2)),'color',[0 0.2 0.6])
% plot(squeeze(nanmean(nanmean(r3,1),2)),'b')
% plot(squeeze(nanmean(nanmean(r_lt0,1),2)),'b')

line([1 25],[0 0],'linestyle',':','color','k')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
axis([1 25 -0.7 0.7])
tp_editplots

subplot(3,2,2); hold on; title('Vertex-wise correlation')
plot(squeeze(r_all1),'color',[0.5 0.2 0])
plot(squeeze(r_all2),'color',[0 0.2 0.6])
% plot(squeeze(r_all3),'color',[0 0.2 1])
% plot(squeeze(r_all4),'color',[0 0.2 0.6])
line([1 25],[0 0],'linestyle',':','color','k')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
axis([1 25 -.6 .6])
tp_editplots
% plot(r_pos_gt0,'r'); plot(abs(r_neg_gt0),'b');
% plot(r_pos_lt0,'r'); plot(abs(r_neg_lt0),'b');

% line([1 25],[0 0],'linestyle',':','color','k')
% set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
% axis([1 25 0 0.60])
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_400grid_splitbydrugeffect_v%d.pdf',v))


%% ALTERED CORR FOR SUBJECT GROUPS

para.nfreq = 1:25; para.alpha = 0.05;
emp1 = pupmod_compute_altered_correlations(cleandat(:,:,SUBJ1,:,:,:),para);
emp2 = pupmod_compute_altered_correlations(cleandat(:,:,SUBJ2,:,:,:),para);

fc1 = squeeze(nanmean(nanmean(cleandat(:,:,SUBJ1,2,2,:),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,SUBJ1,1,2,:),1),2));
fc2 = squeeze(nanmean(nanmean(cleandat(:,:,SUBJ2,2,2,:),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,SUBJ2,1,2,:),1),2));

behav1 = nanmean(behav_cnt(2,SUBJ1,:),3)-nanmean(behav_cnt(1,SUBJ1,:),3);
behav2 = nanmean(behav_cnt(2,SUBJ2,:),3)-nanmean(behav_cnt(1,SUBJ2,:),3);
% behav3 = nanmean(behav_cnt(2,SUBJ3,:),3)-nanmean(behav_cnt(1,SUBJ3,:),3);

para.time = [3000 7000];
pup = pp_loadpupil(para);
pup_cnt = pup(:,:,:,2);

pup_cnt= nanmean(pup_cnt,3);

tmp = pup_cnt(:,2)-pup_cnt(:,1); [i,k]=sort(tmp,'descend');
SUBJ_pup1 = k(1:14); SUBJ_pup2 = k(15:28);

emp_pup1 = pupmod_compute_altered_correlations(cleandat(:,:,SUBJ_pup1,:,:,:),para);
emp_pup2 = pupmod_compute_altered_correlations(cleandat(:,:,SUBJ_pup2,:,:,:),para);

fc_pup1 = squeeze(nanmean(nanmean(cleandat(:,:,SUBJ_pup1,2,2,:),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,SUBJ_pup1,1,2,:),1),2));
fc_pup2 = squeeze(nanmean(nanmean(cleandat(:,:,SUBJ_pup2,2,2,:),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,SUBJ_pup2,1,2,:),1),2));

% tmp = reshape(pup_cnt,[28 3*2]);
% tmp = (tmp-nanmean(tmp,2))./nanstd(tmp,[],2);
% pup_cnt = reshape(tmp,[28 3 2]);

delta_behav1  = nanmean(behav_cnt(2,SUBJ1,:),3)-nanmean(behav_cnt(1,SUBJ1,:),3);
delta_pup1    = pup_cnt(SUBJ1,2)-pup_cnt(SUBJ1,1);
delta_behav2  = nanmean(behav_cnt(2,SUBJ2,:),3)-nanmean(behav_cnt(1,SUBJ2,:),3);
delta_pup2    = pup_cnt(SUBJ2,2)-pup_cnt(SUBJ2,1);
%%
col = cbrewer('seq', 'Reds', 11);

figure; set(gcf,'color','w'); hold on

subplot(4,2,1); hold on
plot(100*emp1.n_p_atx(:,2),'color',col(5,:)); 
plot(100*emp2.n_p_atx(:,2),'color',col(8,:)); 
% plot(100*emp3.n_p_atx(:,2),'color',col(9,:)); 

set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
xlabel('Carrier frequency'); ylabel('Fraction of altered corr. [%]')
axis([1 25 -2 60]); tp_editplots

subplot(4,2,2); hold on
plot(mean(fc1,1),'color',col(5,:)); 
plot(mean(fc2,1),'color',col(8,:)); 
% plot(mean(fc3,1),'color',col(9,:)); 
line([0 25],[0 0],'linestyle',':')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
xlabel('Carrier frequency'); ylabel('\DeltaFC (Atx-Pbo)')
axis([1 25 -0.02 0.04]); tp_editplots
% 

subplot(4,2,3); hold on
bar([1 3],[mean(pup_cnt(SUBJ1,1)) mean(pup_cnt(SUBJ2,1))],0.3,'b')
bar([2 4],[mean(pup_cnt(SUBJ1,2)) mean(pup_cnt(SUBJ2,2))],0.3,'r')
axis([0.5 4.5 7000 8500]); tp_editplots

subplot(4,2,4); hold on
s1 = std((pup_cnt(SUBJ1,1)+pup_cnt(SUBJ1,2))./2)/sqrt(length(SUBJ1));
s2 = std((pup_cnt(SUBJ2,1)+pup_cnt(SUBJ2,2))./2)/sqrt(length(SUBJ2));
m1 = (mean(pup_cnt(SUBJ1,1))+mean(pup_cnt(SUBJ1,2)))./2;
m2 = (mean(pup_cnt(SUBJ2,1))+mean(pup_cnt(SUBJ2,2)))./2;
bar([1],[m1],0.3,'b')
bar([1.5],[m2],0.3,'r')
line([1 1],[m1-s1 m1+s1],'color','k')
line([1.5 1.5],[m2-s2 m2+s2],'color','k')

axis([0.5 2 6500 8500]); tp_editplots


subplot(4,2,7); hold on
plot(100*emp_pup1.n_p_atx(:,2),'color',col(5,:)); 
plot(100*emp_pup2.n_p_atx(:,2),'color',col(8,:)); 
% plot(100*emp3.n_p_atx(:,2),'color',col(9,:)); 

set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
xlabel('Carrier frequency'); ylabel('Fraction of altered corr. [%]')
axis([1 25 -2 60]); tp_editplots


subplot(4,2,8); hold on
plot(mean(fc_pup1,1),'color',col(5,:)); 
plot(mean(fc_pup2,1),'color',col(8,:)); 
% plot(mean(fc3,1),'color',col(9,:)); 
line([0 25],[0 0],'linestyle',':')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
xlabel('Carrier frequency'); ylabel('\DeltaFC (Atx-Pbo)')
axis([1 25 -0.02 0.04]); tp_editplots
% 
% 
% subplot(3,2,3); hold on
% plot(mean(fcbl1,1),'color',col(3,:)); 
% plot(mean(fcbl2,1),'color',col(6,:)); 
% plot(mean(fcbl3,1),'color',col(9,:)); 
% 
% set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128]);
% xlabel('Carrier frequency'); ylabel('\DeltaFC (Atx-Pbo)')
% axis([1 25 -0.02 0.16]); tp_editplots
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fractaltered_splitbygroups_v%d.pdf',v))


figure; set(gcf,'color','w');

subplot(2,2,1); hold on

plot(delta_pup1,delta_behav1,'.','markersize',30,'markerfacecolor','w','markeredgecolor','k'); 
axis square; lsline; tp_editplots; axis([-4000 4000 0 50])
[r1,p1]=corr(delta_pup1(:),delta_behav1(:))
xlabel('\DeltaPupil'); ylabel('\DeltaSR'); title('Group1')
text(-3000,40,sprintf('r = %.3f\np = %.3f',r1,p1),'fontsize',8)

subplot(2,2,2); hold on

plot(delta_pup2,delta_behav2,'.','markersize',30,'markerfacecolor','w','markeredgecolor','k'); 
axis square; lsline; tp_editplots; axis([-4000 4000 -40 10])
[r2,p2]=corr(delta_pup2(:),delta_behav2(:))
xlabel('\DeltaPupil'); ylabel('\DeltaSR'); title('Group2')
text(1000,-30,sprintf('r = %.3f\np = %.3f',r2,p2),'fontsize',8)
% 
% subplot(2,2,3); hold on
% 
% scatter(delta_behav1,mean(fc1(:,10:12),2),100,'.','markerfacecolor','w','markeredgecolor','k'); 
% axis square; lsline; tp_editplots; axis([0 60 -0.04 0.06])
% [r1,p1]=corr(delta_behav1(:),mean(fc1(:,10:12),2))
% xlabel('\DeltaSR'); ylabel('\DeltaFC'); title('Group1')
% text(-3000,40,sprintf('r = %.3f\np = %.3f',r1,p1),'fontsize',8)
% c
% subplot(2,2,4); hold on
% 
% scatter(delta_behav2,mean(fc2(:,10:12),2),100,'.','markerfacecolor','w','markeredgecolor','k'); 
% axis square; lsline; tp_editplots; axis([-40 20 -0.05 0.1])
% [r2,p2]=corr(delta_behav2(:),mean(fc2(:,10:12),2))
% xlabel('\DeltaSR'); ylabel('\DeltaFC'); title('Group2')
% text(1000,-30,sprintf('r = %.3f\np = %.3f',r2,p2),'fontsize',8)


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_400grid_splitbygroups_scatterplot_v%d.pdf',v))

%% PHARMA EFFECT ON SWITCHES - FINE GRAINED



 

%% PERMUTATION TEST
% idx = lab==1
% idx=ones(400,1)==1;
mask    = logical(triu(ones(sum(idx),sum(idx)),1));

nperm = 500;

all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);

behav = behav_cnt(1:2,:,:);

for iperm = 1 : nperm
  iperm
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    behav_perm(1,i,:) = squeeze(behav(idx1(i),i,:));
    behav_perm(2,i,:) = squeeze(behav(idx2(i),i,:));
  end
  
  d_behav_perm   = nanmean(nanmean(behav_perm(2,:,iblock),3),1)-nanmean(nanmean(behav_perm(1,:,iblock),3),1);
  
  [i,k]=sort(d_behav_perm);

  SUBJ1 = k(15:end);

%   SUBJ1 = k(end-sum(d_behav2>0)+1:end);%d_behav2>0;
  SUBJ2 = k(1:14);%d_behav2<0;

  ipharm = 1;

  d_behav_perm    = nanmean(nanmean(behav_perm(ipharm,SUBJ1,iblock),3),1);%-nanmean(nanmean(behav_cnt(1,SUBJ,iblock),3),1);
  d_behav_perm    = permute(repmat(d_behav_perm(:),[1 sum(idx) sum(idx) 25]),[2 3 4 1]);
  d_fc_perm       = squeeze(nanmean(cleandat(idx,idx,SUBJ1,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
  d_fc_perm       = permute(d_fc_perm,[1 2 4 3]);

  % [r2,p2] = tp_corr(d_fc,d_behav,4);
  % [r_all2,p_all2] = tp_corr(squeeze(nanmean(nanmean(d_fc,1),2)),squeeze(nanmean(nanmean(d_behav,1),2)),2);

  % [r1,p1] = tp_corr(d_fc1,d_behav1,4);
  [r_all_perm(:,iperm)] = tp_corr(squeeze(nanmean(nanmean(d_fc_perm,1),2)),squeeze(nanmean(nanmean(d_behav_perm,1),2)),2);

%   for ifoi = 1:25
% 
%     tmp = r2(:,:,ifoi);
%     tmp = tmp((r2(:,:,ifoi)>0));
%     r_pos2(ifoi) = tanh(mean(atanh(tmp)));
% 
%     tmp = r2(:,:,ifoi);
%     tmp = tmp((r2(:,:,ifoi)<0));
%     r_neg2(ifoi) = tanh(mean(atanh(tmp)));
% 
%     tmp = r1(:,:,ifoi);
%     tmp = tmp((r1(:,:,ifoi)>0));
%     r_pos1(ifoi) = tanh(mean(atanh(tmp)));
% 
%     tmp = r1(:,:,ifoi);
%     tmp = tmp((r1(:,:,ifoi)<0));
%     r_neg1(ifoi) = tanh(mean(atanh(tmp)));
%   end

end
%%
p=1-sum(r_all_gt0>r_all_perm(:,1:nperm),2)./nperm;
ci95 = prctile(r_all_perm',95);
ci5 = prctile(r_all_perm',5);

h=p<0.05;
% h=p<fdr1(p,0.05)
% h=(p.*21)<=0.05
figure; set(gcf,'color','w');
subplot(3,2,1); hold on
plot(r_all_gt0,'k')
plot(ci95,'r','linestyle',':')
plot(ci5,'r','linestyle',':')
plot(find(h),r_all1(h),'.','markersize',15)
line([0 26],[0 0],'linestyle',':','color','k')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([1 21 -0.55 0.55]); tp_editplots

p=1-sum(r_all_lt0<r_all_perm(:,1:nperm),2)./nperm;
ci95 = prctile(r_all_perm',95);
ci5 = prctile(r_all_perm',5);

%%

% d_fc       = squeeze(nanmean(cleandat(idx,idx,:,ipharm,cond,:),4));%-squeeze(nanmean(cleandat(idx,idx,SUBJ,1,cond,:),4));
a=nanmean(behav_cnt(2,SUBJ1,:),3)-nanmean(behav_cnt(1,SUBJ1,:),3);
b=squeeze(nanmean(nanmean(nanmean(cleandat(:,:,SUBJ1,2,2,find(h)),6)-nanmean(cleandat(:,:,SUBJ1,1,2,find(h)),6),1),2));

corr(a',b)
%%

%%
  para      = [];
  para.cmap = plasma;
  para.grid = grid;
  para.dd   = 0.75;
  para.clim = [0.1 0.4];
  para.fn   = sprintf('~/pupmod/plots/test.png')
  tp_plot_surface(nanmean(r1(:,:,24),2),para)



%%
clear pp_perm rr_perm r p pos_c*

mask  = logical(tril(ones(400,400),-1));
NPERM = 1000;
ipharm = 1:3;
iblock = 1 :2;
cond   = 2;

if ~exist(sprintf('~/pupmod/proc/conn/pupmod_plots_behav_400grid_permstat_pharm%s.mat',regexprep(num2str(ipharm),' ','')))
  for iperm = 1 : NPERM
    
    fn = sprintf('pupmod_plots_behav_400grid_permstat_ipharm%s_iperm%d',regexprep(num2str(ipharm),' ',''),iperm);
    if tp_parallel(fn,'~/pupmod/proc/conn/',1,0)
      continue
    end
    iperm
    order1 = randperm(28);
    order2 = randperm(28);
    
    d_behav   = nanmean(nanmean(behav_cnt(ipharm,order1,iblock),3),1);
    d_behav   = permute(repmat(d_behav(:),[1 400 400 13]),[2 3 1 4]);
    d_fc      = squeeze(nanmean(nanmean(fc(:,:,order2,ipharm,cond,:,iblock),7),4));
    
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
      cnt = cnt + 1;
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

figure; set(gcf,'color','w');

subplot(2,1,1); hold on

plot(pos_corr,'r');
plot(neg_corr,'b');

tp_editplots;
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
xlabel('Carrier frequency [Hz]'); ylabel('Fraction of sign. correlations [%]')

error('!')


%% LOAD BEHAVIORAL DATA AGAIN


% cleandat = cat(4,cleandat(:,:,:,:,:,:,1),cleandat(:,:,:,:,:,:,2));
% cleandat = (cleandat-nanmean(cleandat,4))./nanstd(cleandat,[],4);
% cleandat = cat(7,cleandat(:,:,:,1:3,:,:),cleandat(:,:,:,4:6,:,:));
% %

para.str_behav  = 'count';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt       = behav;

% 
tmp = reshape(behav_cnt,[3*2 28]);
tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
behav_cnt = reshape(tmp,[3 28 2]);
% 
% para.str_behav  = 'numb_switches';
% behav           = pconn_read_behavioral_data(SUBJLIST,para);
% behav_bttn      = behav;
% behav_bttn      = permute(behav_bttn,[2 1 3]);
% %
% tmp = reshape(behav_bttn,[3*2 28]);
% tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
% behav_bttn = reshape(tmp,[3 28 2]);
% 
% behav_cnt(isnan(behav_cnt))   = behav_bttn(isnan(behav_cnt));
% behav_bttn(isnan(behav_bttn)) = behav_cnt(isnan(behav_bttn));
% 
% behav_cnt = (behav_cnt+behav_bttn)./2;


%% SORT BY BLOCKS (WITHIN CONDITION)
% ---------------------------------
% Identify possible correlations between FC and switch rate by contrasting
% blocks with high switch rate with blocks with low switch rate (within
% subjects)
mask    = logical(triu(ones(400,400),1));
clear s1 s2 r p rr rrr p r d_beh d_fc
ipharm = 1;

SUBJ = 1 : 28
% SUBJ = find(squeeze(nanmean(nanmean(d_fc,1),2))>0);
% if any(any(isnan(mean(behav_cnt,3)'),2))
%   SUBJ(any(isnan(mean(behav_cnt,3)'),2)) = [];
% end

% FOI = [10:12]
 

for FOI = 1:21
%   
  cond = 2; FOI
  s1 = zeros(400,400,length(SUBJ));
  s2 = s1;
    
  [~,idx]=sort(squeeze(nanmean(behav_cnt(ipharm,:,:),1)),2);
  
  for isubj = 1:length(SUBJ)
    s1(:,:,isubj) = (nanmean(nanmean(cleandat(:,:,isubj,ipharm,cond,FOI,idx(isubj,:)==1),6),4));
    s2(:,:,isubj) = (nanmean(nanmean(cleandat(:,:,isubj,ipharm,cond,FOI,idx(isubj,:)==2),6),4));
  end
    
  [h,~,~,s]=ttest(s2,s1,'dim',3);
  
  altered_pos(FOI) = 100*sum((h(mask)>0) & (s.tstat(mask)>0))/sum(mask(:));
  altered_neg(FOI) = 100*sum((h(mask)>0) & (s.tstat(mask)<0))/sum(mask(:));
    
  d_fc=squeeze(nanmean(nanmean(cleandat(:,:,:,2,2,FOI,:),7),6))-squeeze(nanmean(nanmean(cleandat(:,:,:,1,2,FOI,:),7),6));
  d_beh = s2-s1;
  
  mmm(FOI) = nanmean(d_beh(:));
  
  for isubj = 1:length(SUBJ)
    tmp1 = d_fc(:,:,SUBJ(isubj));  tmp1 = tmp1(mask);
    tmp2 = d_beh(:,:,SUBJ(isubj)); tmp2 = tmp2(mask);
    [rr(FOI,isubj)]=corr(tmp1,tmp2,'type','pearson');
%     rmsd(isubj,FOI) = sqrt(mean((nanmean(nanmean(d_fc(:,:,isubj)),3)-nanmean(nanmean(d_beh(:,:,isubj)),3)).^2));

  end
      pattern(:,:,FOI,:)=d_beh;

    [h,~,~,s_fc]=ttest(zeros(400,400,28),d_fc,'dim',3);

%   [rrrrr(FOI), pp(FOI)] = corr(s_fc.tstat(mask),s.tstat(mask));
%   [r(FOI),p(FOI),~,ss]=ttest(zeros(1,size(rr,2)),atan(rr(FOI,:)),'tail','right');
end

%%
figure; set(gcf,'color','w');
subplot(4,4,[2 3]); hold on
plot(nanmean(rr,2),'k')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Correlation coeff.')
line([0 21],[0 0],'linestyle',':','color',[0.7 0.7 0.7])
axis([0 21 -0.1 0.1])

par = nanmean(rr,2);
[~,p]=ttest(zeros(size(rr)),rr,'dim',2,'tail','left')
plot(find(p<0.05),par(p<0.05),'.','markersize',20,'color','k')
tp_editplots


subplot(4,4,[9 10 13 14]); 
h=ttest(zeros(400,400,28),squeeze(pattern(:,:,11,:)),'dim',3);

imagesc(nanmean(pattern(:,:,11,:),4),[-0.02 0.02]); axis off square;
colormap(cmap);title ('High SR - Low SR')
tp_colorbar

d_fc=squeeze(nanmean(nanmean(cleandat(:,:,:,2,2,11,:),7),6))-squeeze(nanmean(nanmean(cleandat(:,:,:,1,2,11,:),7),6));
h=ttest(zeros(size(d_fc)),d_fc,'dim',3)

subplot(4,4,[11 12 15 16]); 
imagesc(nanmean(d_fc,3),[-0.02 0.02]); axis off square;
colormap(cmap);title ('Atx - Pbo')
tp_colorbar


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_behav400_behavior_v%d.pdf',v))

%%
ifoi = 18; 

par = squeeze(nanmean(pattern(:,:,ifoi,:),2));
h = ttest(zeros(size(par)),par,'dim',2);
% for isubj = 1 : 28
para      = [];
para.cmap = cmap;
para.grid = grid;
para.clim = [-0.01 0.01];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_highvslowSR_f%d_pharm%s_v%d.png',ifoi,regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(nanmean(par,2),para)
%
% d_fc=squeeze(nanmean(nanmean(fc(:,:,:,2,2,ifoi,:),7),6))-squeeze(nanmean(nanmean(fc(:,:,:,1,2,FOI,:),7),6));
% d_beh = s2-s1;
% ifoi = 11; 
d_fc=squeeze(nanmean(nanmean(cleandat(:,:,:,2,2,ifoi,:),7),6))-squeeze(nanmean(nanmean(cleandat(:,:,:,1,2,ifoi,:),7),6));

par = squeeze(nanmean(d_fc,2));
h = ttest(zeros(size(par)),par,'dim',2);
% for isubj = 1 : 28
para      = [];
para.cmap = cmap;
para.grid = grid;
para.clim = [-0.01 0.01];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_atxvspbo_f%d_pharm%s_v%d.png',ifoi,regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(nanmean(par,2),para)
% tp_plot_surface(nanmean(nanmean(d_beh(:,:,isubj)),3),para)


% tp_plot_surface(nanmean(nanmean(d_fc(:,:,isubj)),3),para)

% r(isubj)=corr(nanmean(nanmean(d_fc(:,:,isubj)),3)',nanmean(nanmean(d_beh(:,:,isubj)),3)');

% rmsd(isubj) = sqrt(mean((nanmean(nanmean(d_fc(:,:,isubj)),3)-nanmean(nanmean(d_beh(:,:,isubj)),3)).^2));

% end
%% SORT BY BLOCKS (BASED ON DRUG DIFFERENCES)
% ---------------------------------
% Identify possible correlations between FC and switch rate by contrasting
% blocks with high switch rate with blocks with low switch rate (within
% subjects)

clear s1 s2 r p rr rrr p r
% pharm = 1:3;

for FOI = 1:21
  
  cond = 2; FOI
  s1 = zeros(400,400,length(SUBJ));
  s2 = s1;
    
  [~,idx]=sort(squeeze(nanmean(behav_cnt(2,:,:),1))-squeeze(nanmean(behav_cnt(1,:,:),1)),2);
  
  for isubj = 1:28
    s1(:,:,isubj) = s1(:,:,isubj)+(nanmean(nanmean(cleandat(:,:,isubj,ipharm,cond,FOI,idx(isubj,:)==1),6),4))./length(pharm);
    s2(:,:,isubj) = s2(:,:,isubj)+(nanmean(nanmean(cleandat(:,:,isubj,ipharm,cond,FOI,idx(isubj,:)==2),6),4))./length(pharm);
  end
 
  [h,~,~,s]=ttest(s2,s1,'dim',3);
  
  altered_pos(FOI) = sum((h(mask)>0) & (s.tstat(mask)>0));
  altered_neg(FOI) = sum((h(mask)>0) & (s.tstat(mask)<0));
    
  d_fc=squeeze(nanmean(nanmean(cleandat(:,:,SUBJ,2,2,FOI,:),7),6))-squeeze(nanmean(nanmean(cleandat(:,:,SUBJ,1,2,FOI,:),7),6));
  d_beh = s2-s1;
  
  for isubj = 1:size(d_fc,3)
    tmp1 = d_fc(:,:,isubj);  tmp1 = tmp1(mask);
    tmp2 = d_beh(:,:,isubj); tmp2 = tmp2(mask);
    [rr(FOI,isubj)]=corr(tmp1,tmp2,'type','spearman');
  end
  
  [r(FOI),p(FOI),~,ss]=ttest(zeros(1,size(rr,2)),atan(rr(FOI,:)),'tail','right');
  % end
end

%% SORT BY BLOCKS (BASED ON DRUG DIFFERNECES)
% ---------------------------------
% Identify possible correlations between FC and switch rate by contrasting
% blocks with high switch rate with blocks with low switch rate (within
% subjects)

clear s1 s2

ipharm = 2; cond = 2; FOI = 9;

[~,idx]=sort(squeeze(behav_cnt(ipharm,:,:))-squeeze(behav_cnt(1,:,:)),2);

for isubj = 1 : 28
  s1(:,:,isubj,:,:) = fc(:,:,isubj,ipharm,cond,FOI,idx(isubj,1))-fc(:,:,isubj,1,cond,FOI,idx(isubj,1));
  s2(:,:,isubj,:,:) = fc(:,:,isubj,ipharm,cond,FOI,idx(isubj,2))-fc(:,:,isubj,1,cond,FOI,idx(isubj,2));
end

[h,p,~,s]=ttest(squeeze(nanmean(nanmean(s2,1),5)),squeeze(nanmean(nanmean(s1,1),5)),'dim',2);

para      = [];
para.cmap = cmap;
para.grid = grid;
para.clim = [-2 2];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_sortbyblocks_f%d_pharm%s_v%d.png',ifoi,regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(squeeze(s.tstat),para)

para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_sortbyblocks_masked_f%d_pharm%s_v%d.png',ifoi,regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(squeeze(s.tstat).*(p<0.05),para)


%% SORT BY BLOCKS: SPECTRUM
clear s1 s2
ipharm = 2; cond = 2; FOI = 22:22;

[~,idx]=sort(squeeze(behav_cnt(ipharm,:,:)),2);

for isubj = 1 : 28
  s1(:,:,isubj,:,:) = fc(:,:,isubj,:,cond,:,idx(isubj,1));
  s2(:,:,isubj,:,:) = fc(:,:,isubj,:,cond,:,idx(isubj,2));
end

[h,p,~,s]=ttest(squeeze(nanmean(s2(:,:,:,ipharm,:),1)),squeeze(nanmean(s1(:,:,:,ipharm,:),1)),'dim',2);

pos = sum((squeeze(s.tstat)>0)&(squeeze(h)>0));
neg = sum((squeeze(s.tstat)<0)&(squeeze(h)>0));

%% SORT BY BLOCKS: PERMUTATION TEST

ipharm = 2; cond = 2; FOI = 22:22;

tmp = squeeze(fc(:,:,:,ipharm,cond,:,:));

for iperm = 1 : 1500
  
  iperm
  idx(:,1) = randi(2,[28 1]);
  idx(:,2) = 3-idx(:,1);
  
  for isubj = 1 : 28
    s1(:,:,isubj,:) = tmp(:,:,isubj,:,idx(isubj,1));
    s2(:,:,isubj,:) = tmp(:,:,isubj,:,idx(isubj,2));
  end
  
  [h,p,~,s]=ttest(squeeze(nanmean(s2,1)),squeeze(nanmean(s1,1)),'dim',2);
  
  pos_perm(:,iperm) = sum((squeeze(s.tstat)>0)&(squeeze(h)>0));
  neg_perm(:,iperm) = sum((squeeze(s.tstat)<0)&(squeeze(h)>0));
  
end


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

%%

for ifoi = 1 : 25
  
  [h,~,~,s]=ttest(nanmean(fc(:,:,:,2,2,ifoi,:),7),nanmean(fc(:,:,:,1,2,ifoi,:),7),'dim',3);
  
  pos(ifoi) = sum((h(mask)>0) & (s.tstat(mask)>0));
  neg(ifoi) = sum((h(mask)>0) & (s.tstat(mask)<0));
  
end

