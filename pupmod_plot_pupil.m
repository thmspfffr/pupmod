%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

%% INFORMATION
% First plot everything in VTPM atlas
% See below for whole cortex

clear

v = 20;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/
addpath ~/pp/matlab
outdir = '~/pupmod/proc/conn/';

cleandat = pupmod_loadpowcorr(v,0);

if v == 20
  tmp = tp_create_grid('vtpm');
  idx = 1:46;
  idx = ~ismember(idx,[21 22 23 44 45 46]);
  reg = tmp.tissuelabel_4mm(idx);
end

%% LOAD BEHAVIOR & PUPIL

para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt = behav;

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);

para.time = [3000 7000];
pup = pp_loadpupil(para);

%% BASLINE: correlate placebo (pupil) with placebo (FC)
% AVERAGE OVER BLOCKS

clear r_cnt p_cnt

ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean   = squeeze(nanmean(cleandat(idx,idx,:,:,2,ifoi,:),7));
mask      = logical(tril(ones(40,40),-1));

d_pup     = nanmean(pup(:,1,:,2),3);
d_fc      = pc_mean(:,:,:,1);

% identify and ignore nans
nan_idx_cnt   = ~isnan(d_pup);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_full_v%d.pdf',v))

clear r_cnt p_cnt r_bttn p_bttn

% ------------------------
% LEFT/RIGHT COLLAPSED
% ------------------------
pc_mean   = (squeeze(nanmean(cleandat(1:20,1:20,:,:,2,ifoi,:),7))+squeeze(nanmean(cleandat(24:43,24:43,:,:,2,ifoi,:),7)))./2;
pc_blocks = (squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:))+squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:)))./2;
d_fc      = pc_mean(:,:,:,1);

% compute correlation
for i = 1 :20
  for j = 1 : 20
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_lr_v%d.pdf',v))

clear r_cnt p_cnt r_bttn p_bttn

% ------------------------
% DR MURPHY REGIONS
% ------------------------

labels = {'V1';'V2-4';'V3ab';'IPS01';'IPS23';'LO';'VO';'PHC';'MT'};
region{1} = [1 2]; region{2} = [3 4 5 6 7]; region{3} = [15 16]; region{4} = [17 18];
region{5} = [19 20]; region{6} = [13 14]; region{7} = [8 9]; region{8} = [10 11]; region{9} = [12];

cleandat_lr = squeeze((cleandat(1:20,1:20,:,:,2,6,:)+cleandat(24:43,24:43,:,:,2,6,:))./2); 
clear tmp fc

% collapse across regions defined in 'labels'
for i = 1 : size(region,2)
  for j = 1 : size(region,2)
    if i == j; fc(i,j,:,:,:) = nan(28,3,2); continue; end
    clear tmp
    for ii = 1 : length(region{i})  
      tmp(ii,:,:,:) = squeeze(mean(cleandat_lr(region{i}(ii),region{j},:,:,:),2));
    end
    fc(i,j,:,:,:) = squeeze(mean(tmp,1));    
  end
end

clear r rb p pb
pc_mean   = nanmean(fc,5);
d_fc      = pc_mean(:,:,:,1);

for i = 1 :9
  for j = 1 : 9
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j),p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_coarse_v%d.pdf',v))

%% DO THE SAME SEPARATELY FOR THE TWO BLOCKS
if v==20
mask = logical(tril(ones(40,40),-1));
else
  mask = logical(tril(ones(400,400),-1));
end
  
ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean1   = squeeze(cleandat(idx,idx,:,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(idx,idx,:,:,2,ifoi,2));

d_pup1     = pup(:,1,1,2);
d_pup2     = pup(:,1,2,2);

d_fc1      = pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_pup1); 
nan_idx_cnt2   = ~isnan(d_pup2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1));
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2)')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_full_blocks_v%d.pdf',v))

clear r_cnt1 p_cnt1 r_cnt2 p_cnt2

% ------------------------
% LEFT/RIGHT COLLAPSED
% ------------------------
pc_blocks  = (squeeze(cleandat(1:20,1:20,:,:,2,ifoi,:))+squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:)))./2;

pc_mean1   = squeeze(pc_blocks(:,:,:,:,1));
pc_mean2   = squeeze(pc_blocks(:,:,:,:,2));

d_fc1      = pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,1);

% compute correlation
for i = 1 : 20
  for j = 1 : 20
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1));
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2)')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_lr_blocks_v%d.pdf',v))

clear r_cnt1 p_cnt1 r_cnt2 p_cnt2

% ------------------------
% DR MURPHY REGIONS
% ------------------------

labels = {'V1';'V2-4';'V3ab';'IPS01';'IPS23';'LO';'VO';'PHC';'MT'};
region{1} = [1 2]; region{2} = [3 4 5 6 7]; region{3} = [15 16]; region{4} = [17 18];
region{5} = [19 20]; region{6} = [13 14]; region{7} = [8 9]; region{8} = [10 11]; region{9} = [12];

cleandat_lr = squeeze((cleandat(1:20,1:20,:,:,2,6,:)+cleandat(24:43,24:43,:,:,2,6,:))./2); 
clear tmp fc

% collapse across regions defined in 'labels'
for i = 1 : size(region,2)
  for j = 1 : size(region,2)
    if i == j; fc(i,j,:,:,:) = nan(28,3,2); continue; end
    clear tmp
    for ii = 1 : length(region{i})  
      tmp(ii,:,:,:) = squeeze(mean(cleandat_lr(region{i}(ii),region{j},:,:,:),2));
    end
    fc(i,j,:,:,:) = squeeze(mean(tmp,1));    
  end
end

pc_mean1   = squeeze(fc(:,:,:,:,1));
pc_mean2   = squeeze(fc(:,:,:,:,2));

d_fc1      = pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,1);

for i = 1 :9
  for j = 1 : 9
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j),p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1));
    [r_cnt2(i,j),p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square; title('Counting: Block #2)')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_coarse_blocks_v%d.pdf',v))


%% CORRELATE FC difference WITH PUPIL difference: data averaged over blocks

clear r_cnt p_cnt

ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean   = squeeze(nanmean(cleandat(idx,idx,:,:,2,ifoi,:),7));
mask      = logical(tril(ones(40,40),-1));

d_pup     = nanmean(pup(:,2,:,2),3)-nanmean(pup(:,1,:,2),3);
d_fc      = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

% identify and ignore nans
nan_idx_cnt   = ~isnan(d_pup);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_full_v%d.pdf',v))

clear r_cnt p_cnt r_bttn p_bttn

% ------------------------
% LEFT/RIGHT COLLAPSED
% ------------------------
pc_mean   = (squeeze(nanmean(cleandat(1:20,1:20,:,:,2,ifoi,:),7))+squeeze(nanmean(cleandat(24:43,24:43,:,:,2,ifoi,:),7)))./2;
pc_blocks = (squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:))+squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:)))./2;
d_fc      = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

% compute correlation
for i = 1 :20
  for j = 1 : 20
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_lr_v%d.pdf',v))

clear r_cnt p_cnt r_bttn p_bttn

% ------------------------
% DR MURPHY REGIONS
% ------------------------

labels = {'V1';'V2-4';'V3ab';'IPS01';'IPS23';'LO';'VO';'PHC';'MT'};
region{1} = [1 2]; region{2} = [3 4 5 6 7]; region{3} = [15 16]; region{4} = [17 18];
region{5} = [19 20]; region{6} = [13 14]; region{7} = [8 9]; region{8} = [10 11]; region{9} = [12];

cleandat_lr = squeeze((cleandat(1:20,1:20,:,:,2,6,:)+cleandat(24:43,24:43,:,:,2,6,:))./2); 
clear tmp fc

% collapse across regions defined in 'labels'
for i = 1 : size(region,2)
  for j = 1 : size(region,2)
    if i == j; fc(i,j,:,:,:) = nan(28,3,2); continue; end
    clear tmp
    for ii = 1 : length(region{i})  
      tmp(ii,:,:,:) = squeeze(mean(cleandat_lr(region{i}(ii),region{j},:,:,:),2));
    end
    fc(i,j,:,:,:) = squeeze(mean(tmp,1));    
  end
end

clear r rb p pb
pc_mean   = nanmean(fc,5);
d_fc      = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

for i = 1 :9
  for j = 1 : 9
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j),p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_coarse_v%d.pdf',v))

%% DO THE SAME SEPARATELY FOR THE TWO BLOCKS

corr_type = 'spearman';

mask = logical(tril(ones(40,40),-1));
ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean1   = squeeze(cleandat(idx,idx,:,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(idx,idx,:,:,2,ifoi,2));

d_pup1     = pup(:,2,1,2)-pup(:,1,1,2);
d_pup2     = pup(:,2,2,2)-pup(:,1,2,2);

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_pup1); 
nan_idx_cnt2   = ~isnan(d_pup2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1),'type',corr_type);
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1),'type',corr_type);
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2)')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_full_blocks_v%d.pdf',v))

clear r_cnt1 p_cnt1 r_cnt2 p_cnt2

% ------------------------
% LEFT/RIGHT COLLAPSED
% ------------------------
pc_blocks  = (squeeze(cleandat(1:20,1:20,:,:,2,ifoi,:))+squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:)))./2;

pc_mean1   = squeeze(pc_blocks(:,:,:,:,1));
pc_mean2   = squeeze(pc_blocks(:,:,:,:,2));

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

% compute correlation
for i = 1 : 20
  for j = 1 : 20
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1),'type',corr_type);
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1),'type',corr_type);
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2)')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_lr_blocks_v%d.pdf',v))

clear r_cnt1 p_cnt1 r_cnt2 p_cnt2

% ------------------------
% DR MURPHY REGIONS
% ------------------------

labels = {'V1';'V2-4';'V3ab';'IPS01';'IPS23';'LO';'VO';'PHC';'MT'};
region{1} = [1 2]; region{2} = [3 4 5 6 7]; region{3} = [15 16]; region{4} = [17 18];
region{5} = [19 20]; region{6} = [13 14]; region{7} = [8 9]; region{8} = [10 11]; region{9} = [12];

cleandat_lr = squeeze((cleandat(1:20,1:20,:,:,2,6,:)+cleandat(24:43,24:43,:,:,2,6,:))./2); 
clear tmp fc

% collapse across regions defined in 'labels'
for i = 1 : size(region,2)
  for j = 1 : size(region,2)
    if i == j; fc(i,j,:,:,:) = nan(28,3,2); continue; end
    clear tmp
    for ii = 1 : length(region{i})  
      tmp(ii,:,:,:) = squeeze(mean(cleandat_lr(region{i}(ii),region{j},:,:,:),2));
    end
    fc(i,j,:,:,:) = squeeze(mean(tmp,1));    
  end
end

pc_mean1   = squeeze(fc(:,:,:,:,1));
pc_mean2   = squeeze(fc(:,:,:,:,2));

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

for i = 1 :9
  for j = 1 : 9
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j),p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_pup1(nan_idx1),[],1),'type',corr_type);
    [r_cnt2(i,j),p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_pup2(nan_idx2),[],1),'type',corr_type);
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1)')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square; title('Counting: Block #2)')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_vtpm_coarse_blocks_v%d.pdf',v))

%% SCATTER PLOTS: Change in FC (avg across all VTPM) vs. change in pupil

mask = logical(tril(ones(40,40),-1));
ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean1   = squeeze(cleandat(idx,idx,:,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(idx,idx,:,:,2,ifoi,2));

d_pup1     = pup(:,2,1,2)-pup(:,1,1,2);
d_pup2     = pup(:,2,2,2)-pup(:,1,2,2);

d_pup1_prc = 100*(pup(:,2,1,2)-pup(:,1,1,2))./pup(:,1,1,2);
d_pup2_prc = 100*(pup(:,2,2,2)-pup(:,1,2,2))./pup(:,1,2,2);

d_fc1      = nansum(reshape(pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1),[size(pc_mean1,1)*size(pc_mean1,1) 28]))./(size(pc_mean1,1)*size(pc_mean1,1));
d_fc2      = nansum(reshape(pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1),[size(pc_mean1,1)*size(pc_mean1,1) 28]))./(size(pc_mean1,1)*size(pc_mean1,1));

d_fc1_prc  = 100*nansum(reshape((pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1))./pc_mean1(:,:,:,1),[size(pc_mean1,1)*size(pc_mean1,1) 28]))./(size(pc_mean1,1)*size(pc_mean1,1));
d_fc2_prc  = 100*nansum(reshape((pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1))./pc_mean2(:,:,:,1),[size(pc_mean1,1)*size(pc_mean1,1) 28]))./(size(pc_mean1,1)*size(pc_mean1,1));

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_pup1); 
nan_idx_cnt2   = ~isnan(d_pup2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

figure; set(gcf,'color','w');

subplot(2,2,1);
scatter(d_pup1(nan_idx1),d_fc1(nan_idx1))
[r,p]=corr(reshape(d_pup1(nan_idx1),[],1),reshape(d_fc1(nan_idx1),[],1));
axis square; tp_editplots; xlabel('Change in pupil diameter');  ylabel('Change in mean FC')
axis([ -4000 4000 -0.2 0.2]); lsline
text(-3000,0.13,sprintf('r = %.3f\np = %.3f',r,p),'fontsize',7)

subplot(2,2,2);
scatter(d_pup2(nan_idx2),d_fc2(nan_idx2))
[r,p]=corr(reshape(d_pup2(nan_idx2),[],1),reshape(d_fc2(nan_idx2),[],1))
axis square; tp_editplots;  xlabel('Change in pupil diameter');  ylabel('Change in mean FC')
axis([ -4000 4000 -0.2 0.2]); lsline
text(-3000,0.13,sprintf('r = %.3f\np = %.3f',r,p),'fontsize',7)

subplot(2,2,3);
scatter(d_pup1_prc(nan_idx1),d_fc1_prc(nan_idx1))
[r,p]=corr(reshape(d_pup1_prc(nan_idx1),[],1),reshape(d_fc1_prc(nan_idx1),[],1));
axis square; tp_editplots; xlabel('Change in pupil diameter [%]');  ylabel('Change in mean FC [%]')
axis([-100 100 -200 200 ]); lsline
text(-80,-100,sprintf('r = %.3f\np = %.3f',r,p),'fontsize',7)

subplot(2,2,4);
scatter(d_pup2_prc(nan_idx2),d_fc2_prc(nan_idx2))
[r,p]=corr(reshape(d_pup2_prc(nan_idx2),[],1),reshape(d_fc2_prc(nan_idx2),[],1))
axis square; tp_editplots;  xlabel('Change in pupil diameter [%]');  ylabel('Change in mean FC [%]')
axis([-100 100 -200 200]); lsline
text(-80,-100,sprintf('r = %.3f\np = %.3f',r,p),'fontsize',7)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_pupil_scatter_vtpm_v%d.pdf',v))

%% JUST PUPIL

figure; set(gcf,'color','w');

m1 = nanmean(pup(:,[1 2],1,2));
m2 = nanmean(pup(:,[1 2],2,2));
s1 = sqrt(nansum((pup(:,[1 2],1,2)-nanmean(pup(:,[1 2],1,2))).^2)/27)/sqrt(28);
s2 = sqrt(nansum((pup(:,[1 2],2,2)-nanmean(pup(:,[1 2],2,2))).^2)/27)/sqrt(28);

subplot(2,2,1);
bar(m1);  title('Block #1')
line([1 1],[m1(1)-s1(1) m1(1)+s1(1)],'linestyle','-','color','k')
line([2 2],[m1(2)-s1(2) m1(2)+s1(2)],'linestyle','-','color','k')
set(gca,'xTick',[1 2],'xTickLabels',{'Pbo';'Atx'},'ticklabelinterpreter','none');xtickangle(90)
tp_editplots; ylabel('Pupil diameter'); axis([0 3 6000 9000])
subplot(2,2,2); 
bar(m2);  title('Block #2')
line([1 1],[m2(1)-s2(1) m2(1)+s2(1)],'linestyle','-','color','k')
line([2 2],[m2(2)-s2(2) m2(2)+s2(2)],'linestyle','-','color','k')
set(gca,'xTick',[1 2],'xTickLabels',{'Pbo';'Atx'},'ticklabelinterpreter','none');xtickangle(90)
tp_editplots; ylabel('Pupil diameter'); axis([0 3 6000 9000])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_bar_v%d.pdf',v))


%% WHOLE CORTEX FOR PLACEBO

clear cleandat

v = 12;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
outdir = '~/pupmod/proc/conn/';
addpath /home/gnolte/meg_toolbox/toolbox_nightly/

cleandat = pupmod_loadpowcorr(v,0);
load sa_meg_template;

fprintf('Loading grid...\n')
grid  = select_chans(sa_meg_template.grid_cortex3000,400); 
fprintf('Loading grid... Done\n')

%% PLOT FC MATRIX FOR DIFFERENT CONDITIONS

clear r_cnt p_cnt

ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
mask      = logical(tril(ones(400,400),-1));

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean   = squeeze(nanmean(cleandat(:,:,:,:,2,ifoi,:),7));
d_pup     = nanmean(pup(:,1,:,2),3);
d_fc      = pc_mean(:,:,:,1);

% identify and ignore nans
nan_idx_cnt   = ~isnan(d_pup);

% compute correlation
for i = 1 :400
  i
  for j = 1 : 400
    if i == j; r_cnt(i,j)=nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_pup(nan_idx_cnt),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_pupil_baseline_vtpm_full_v%d.pdf',v))

clear r_cnt p_cnt r_bttn p_bttn




FOI = 6:7

if v == 12
  for ifoi = FOI

    figure; set(gcf,'color','w')
    cond = [1 2 3];
    [~,front_to_back] = sort(grid(:,2),'descend');
    left  = find(grid(:,1)<0);
    right = find(grid(:,1)>0);
    for i = 1 : 3
      fc = nanmean(nanmean(cleandat(:,:,:,cond(i),1,ifoi,:),7),3);

      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));

      fc_rest = [fc1 fc2; fc3 fc4];

      fc = nanmean(nanmean(cleandat(:,:,:,cond(i),2,ifoi,:),7),3);

      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));

      fc_task = [fc1 fc2; fc3 fc4];

      fc = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(1,3,i); imagesc(fc,[0.02 0.1]); axis square off
      colormap(inferno)

    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_f%d_v%d.pdf',ifoi,v))
  end

end