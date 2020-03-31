%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 3;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/
addpath ~/pp/matlab
outdir = '~/pupmod/proc/conn/';

cleandat = pupmod_loadpowcorr(v,SUBJLIST,0);

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

%% PLOT BEHAVIOR

allpup_cnt  = pup(:,:,:,2);
allpup_bttn = pup(:,:,:,3);

i_pharm = 2;
corr_type = 'spearman';

figure; set(gcf,'color','w'); hold on

for iblock = 1 : 2

  prc_pup_cnt(:,iblock)  = 100*(allpup_cnt(:,i_pharm,iblock)-allpup_cnt(:,1,iblock))./allpup_cnt(:,1,iblock);
  prc_pup_btn(:,iblock)  = 100*(allpup_bttn(:,i_pharm,iblock)-allpup_bttn(:,1,iblock))./allpup_bttn(:,1,iblock);
  prc_beh_cnt(:,iblock) = 100*(behav_cnt(i_pharm,:,iblock)-behav_cnt(1,:,iblock))./behav_cnt(1,:,iblock);

  nanidx_cnt3 = abs(prc_beh_cnt(:,iblock))>100;
  nanidx_cnt1 = isnan(prc_beh_cnt(:,iblock)); nanidx_cnt2 = isnan(prc_pup_cnt(:,iblock));
  nanidx_cnt  = ~(nanidx_cnt1|nanidx_cnt2);

  % PLOT
  % -----------
  subplot(2,2,iblock)
  scatter(prc_beh_cnt(nanidx_cnt),prc_pup_cnt(nanidx_cnt),50,'markerfacecolor','k','markeredgecolor','w');
  [r,p]=corr(prc_beh_cnt(nanidx_cnt,iblock),prc_pup_cnt(nanidx_cnt,iblock),'type',corr_type)
  axis([-50 220 -50 120]); tp_editplots
  text(-38,100,sprintf('r = %.3f | p = %.3f',r,p))
  title(sprintf('Counting Block%d',iblock)); lsline
  xlabel('\DeltaBehavior [%]'); ylabel('\DeltaPupil [%]'); tp_editplots
  % -----------

  prc_beh_bttn(:,iblock) = 100*(behav_bttn(i_pharm,:,iblock)-behav_bttn(1,:,iblock))./behav_bttn(1,:,iblock);
  prc_pup_btn(prc_beh_bttn>300,iblock)  = nan;
  prc_beh_bttn(prc_beh_bttn>300,iblock) = nan;

  nanidx_cnt3 = abs(prc_beh_bttn(:,iblock))>100;

  nanidx_cnt1 = isnan(prc_beh_bttn(:,iblock)); nanidx_cnt2 = isnan(prc_pup_btn(:,iblock));
  nanidx_cnt = ~(nanidx_cnt1|nanidx_cnt2);

  % PLOT
  % -----------
  subplot(2,2,iblock+2)
  scatter(prc_beh_bttn(nanidx_cnt,iblock),prc_pup_btn(nanidx_cnt,iblock),50,'markerfacecolor','k','markeredgecolor','w');
  [r,p]=corr(prc_beh_bttn(nanidx_cnt,iblock),prc_pup_btn(nanidx_cnt,iblock),'type',corr_type)
  axis([-50 250 -50 120]); tp_editplots
  text(-38,100,sprintf('r = %.3f | p = %.3f',r,p))
  title(sprintf('Pressing Block%d',iblock)); lsline
  xlabel('\DeltaBehavior [%]'); ylabel('\DeltaPupil [%]');tp_editplots
  % -----------
end
print(gcf,'-dpdf',sprintf('~/pp/plots/pupil_behav_indivblocks.pdf'))

% POOLED PLOTS
% --------------------------
aa = prc_beh_bttn(:); 
ab = prc_pup_btn(:); 

nanidx = isnan(aa)|isnan(ab);
aa(nanidx) = []; ab(nanidx) = [];

ba = prc_beh_cnt(:);
bb = prc_pup_cnt(:);
nanidx = isnan(ba)|isnan(bb);
ba(nanidx) = []; bb(nanidx) = [];

figure; set(gcf,'color','w'); hold on
scatter([aa; ba],[ab; bb],50,'markerfacecolor','k','markeredgecolor','w');
lsline; axis([-150 250 -80 100]);
line([0 0],[-80 100],'color',[0.8 0.8 0.8],'linestyle',':')
line([-150 250],[0 0],'color',[0.8 0.8 0.8],'linestyle',':')
[r,p]=corr([aa; ba],[ab; bb],'type','spearman')
text(-38,100,sprintf('r = %.3f | p = %.3f',r,p))
xlabel('\DeltaBehavior [%]'); ylabel('\DeltaPupil [%]')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pupil_behav_pooled.pdf'))


%% CORRELATE FC WITH BEHAVIOR: data averaged over blocks
% FINISH THIS, CURRENTLY THROWS ERROR

ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean   = squeeze(cleandat(idx,idx,:,:,2,ifoi,:));
pc_mean   = zscore(reshape(pc_mean,[40 40 28 6]),0,4);
pc_mean   = reshape(pc_mean,[40 40 28 3 2]);
pc_mean   = nanmean(pc_mean,5);
mask      = logical(tril(ones(40,40),-1));

a=reshape(permute(behav_cnt,[2 1 3]),[28 6]);
behav_z = permute(reshape((a - nanmean(a)) ./ nanstd(a),[28 3 2]),[2 1 3]);

% d_behav       = nanmean(behav_cnt(2,:,:)-behav_cnt(1,:,:),3);
% d_behav_bttn  = nanmean(behav_bttn(2,:,:)-behav_bttn(1,:,:),3);
% d_fc          = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

d_behav = nanmean(behav_z(2,:,:)-behav_z(1,:,:),3);
d_fc          = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

% identify and ignore nans
nan_idx_cnt   = ~isnan(d_behav);
nan_idx_fc    = ~isnan(squeeze(d_fc(1,2,:)));
nan_idx_cnt = nan_idx_cnt(:)&nan_idx_fc(:);
% nan_idx_bttn  = ~isnan(d_behav_bttn);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt(i,j)=nan; r_bttn(i,j) = nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_behav(nan_idx_cnt),[],1));
%     [r_bttn(i,j), p_bttn(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_bttn)),reshape(d_behav_bttn(nan_idx_bttn),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w counting')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

% subplot(2,2,3); imagesc(r_bttn,[-0.3 0.3]); axis square; title('Correlation w pressing')
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
% tp_editplots; set(gca,'FontSize',5)
% subplot(2,2,4); imagesc(r_bttn.*(p_bttn<0.05),[-0.3 0.3]); axis square;
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
% colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_full_v%d.pdf',v))
%%
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
    if i == j; r_cnt(i,j)=nan; r_bttn(i,j) = nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_behav(nan_idx_cnt),[],1));
%     [r_bttn(i,j), p_bttn(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_bttn)),reshape(d_behav_bttn(nan_idx_bttn),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w counting')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

% subplot(2,2,3); imagesc(r_bttn,[-0.3 0.3]); axis square; title('Correlation w pressing')
% set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
% tp_editplots; set(gca,'FontSize',5)
% subplot(2,2,4); imagesc(r_bttn.*(p_bttn<0.05),[-0.3 0.3]); axis square;
% set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
% colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_lr_v%d.pdf',v))

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
    if i == j; r_cnt(i,j)=nan; r_bttn(i,j) = nan; continue; end
    
    [r_cnt(i,j),p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_cnt)),reshape(d_behav(nan_idx_cnt),[],1));
%     [r_bttn(i,j),p_bttn(i,j)] = corr(squeeze(d_fc(i,j,nan_idx_bttn)),reshape(d_behav_bttn(nan_idx_bttn),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Correlation w counting')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_coarse_v%d.pdf',v))

%% DO THE SAME SEPARATELY FOR THE TWO BLOCKS
SUBJ = 1:28; %SUBJ(13)=[];

mask = logical(tril(ones(40,40),-1));
ifoi = 6; cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

% START WITH FULL VTPM ATLAS
% ------------------------
pc_mean1   = squeeze(cleandat(idx,idx,SUBJ,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(idx,idx,SUBJ,:,2,ifoi,2));

d_behav1   = behav_cnt(2,SUBJ,1)-behav_cnt(1,SUBJ,1);
d_behav2   = behav_cnt(2,SUBJ,2)-behav_cnt(1,SUBJ,2);

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_behav1); 
nan_idx_cnt2   = ~isnan(d_behav2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

% compute correlation
for i = 1 :40
  for j = 1 : 40
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_behav1(nan_idx1),[],1));
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_behav2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_full_blocks_v%d.pdf',v))

clear r_cnt1 p_cnt1 r_cnt2 p_cnt2

% ------------------------
% LEFT/RIGHT COLLAPSED
% ------------------------
pc_blocks  = (squeeze(cleandat(1:20,1:20,:,:,2,ifoi,:))+squeeze(cleandat(24:43,24:43,:,:,2,ifoi,:)))./2;

pc_mean1   = squeeze(pc_blocks(:,:,:,:,1));
pc_mean2   = squeeze(pc_blocks(:,:,:,:,2));

d_behav1   = behav_cnt(2,:,1)-behav_cnt(1,:,1);
d_behav2   = behav_cnt(2,:,2)-behav_cnt(1,:,2);

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_behav1); 
nan_idx_cnt2   = ~isnan(d_behav2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

% compute correlation
for i = 1 : 20
  for j = 1 : 20
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_behav1(nan_idx1),[],1));
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_behav2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title('Counting: Block #2')
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:20,'xTickLabels',reg(1:1:20),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:20,'yTickLabels',reg(1:1:20),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_lr_blocks_v%d.pdf',v))

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

d_behav1   = behav_cnt(2,:,1)-behav_cnt(1,:,1);
d_behav2   = behav_cnt(2,:,2)-behav_cnt(1,:,2);

d_fc1      = pc_mean1(:,:,:,2)-pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,2)-pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_behav1); 
nan_idx_cnt2   = ~isnan(d_behav2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

for i = 1 :9
  for j = 1 : 9
    if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt1(i,j),p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_behav1(nan_idx1),[],1));
    [r_cnt2(i,j),p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_behav2(nan_idx2),[],1));
  end
end

figure; set(gcf,'color','w');

subplot(2,2,1); imagesc(r_cnt1,[-0.3 0.3]); axis square; title('Counting: Block #1')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,2); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

subplot(2,2,3); imagesc(r_cnt2,[-0.3 0.3]); axis square; title('Counting: Block #2')
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(2,2,4); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:1:9,'xTickLabels',labels(1:1:9),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:9,'yTickLabels',labels(1:1:9),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_vtpm_coarse_blocks_v%d.pdf',v))

%% CORRELATE PLACEBO BEHAVIOR WITH PLACEBO POWCORR
% START WITH FULL VTPM ATLAS
% ------------------------
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ = 1:28; %SUBJ(13)=[];

for ifoi = 1:13

clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2

pc_mean1   = squeeze(cleandat(idx,idx,SUBJ,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(idx,idx,SUBJ,:,2,ifoi,2));
% pc_mean1   = squeeze(fc(:,:,SUBJ,2,2,ifoi,1))-squeeze(fc(:,:,SUBJ,1,2,ifoi,1));
% pc_mean2   = squeeze(fc(:,:,SUBJ,2,2,ifoi,2))-squeeze(fc(:,:,SUBJ,1,1,ifoi,2));

d_behav    = nanmean(behav_cnt(1,SUBJ,:),3);
d_behav1   = behav_cnt(1,SUBJ,1);%behav_cnt(2,SUBJ,1)-behav_cnt(1,SUBJ,1);
d_behav2   = behav_cnt(2,SUBJ,2);%behav_cnt(2,SUBJ,2)-behav_cnt(1,SUBJ,2);

d_fc       = squeeze(nanmean(cleandat(idx,idx,SUBJ,1,2,ifoi,:),7));
% d_fc       = squeeze(nanmean(fc(:,:,SUBJ,2,2,ifoi,:),7))-squeeze(nanmean(fc(:,:,SUBJ,1,2,ifoi,:),7));

d_fc1      = pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_behav1); 
nan_idx_cnt2   = ~isnan(d_behav2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

% compute correlation
for i = 1 : size(d_fc1,1)
  for j = 1 : size(d_fc1,1)
%     if i == j; r_cnt1(i,j)=nan; r_cnt2(i,j) = nan; continue; end
    [r_cnt(i,j), p_cnt(i,j)] = corr(squeeze(d_fc(i,j,nan_idx1)),reshape(d_behav(nan_idx1),[],1));
    [r_cnt1(i,j), p_cnt1(i,j)] = corr(squeeze(d_fc1(i,j,nan_idx1)),reshape(d_behav1(nan_idx1),[],1));
    [r_cnt2(i,j), p_cnt2(i,j)] = corr(squeeze(d_fc2(i,j,nan_idx2)),reshape(d_behav2(nan_idx2),[],1));
  end
end

mask      = logical(tril(ones(size(d_fc1,1),size(d_fc1,1)),-1));
r12(ifoi) = corr(r_cnt1(mask),r_cnt2(mask))

figure; set(gcf,'color','w');

subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Counting: Average')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %d',ifoi))
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %d)',ifoi))
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5);
subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_placebo_vtpm_blocks_f%d_v%d.pdf',ifoi,v))
end

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

[~,front_to_back] = sort(grid(:,2),'descend');
left  = find(grid(:,1)<0);
right = find(grid(:,1)>0);

%% PLOT FC MATRIX FOR DIFFERENT CONDITIONS
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
SUBJ = 1:28; %SUBJ(13)=[];

for ifoi = 1:13

clear r_cnt p_cnt r_cnt1 p_cnt1 r_cnt2 p_cnt2

pc_mean1   = squeeze(cleandat(:,:,SUBJ,:,2,ifoi,1));
pc_mean2   = squeeze(cleandat(:,:,SUBJ,:,2,ifoi,2));
% pc_mean1   = squeeze(fc(:,:,SUBJ,2,2,ifoi,1))-squeeze(fc(:,:,SUBJ,1,2,ifoi,1));
% pc_mean2   = squeeze(fc(:,:,SUBJ,2,2,ifoi,2))-squeeze(fc(:,:,SUBJ,1,1,ifoi,2));

d_behav    = nanmean(behav_cnt(1,SUBJ,:),3);
d_behav1   = behav_cnt(1,SUBJ,1);%behav_cnt(2,SUBJ,1)-behav_cnt(1,SUBJ,1);
d_behav2   = behav_cnt(2,SUBJ,2);%behav_cnt(2,SUBJ,2)-behav_cnt(1,SUBJ,2);

d_fc       = squeeze(nanmean(cleandat(:,:,SUBJ,1,2,ifoi,:),7));

d_fc1      = pc_mean1(:,:,:,1);
d_fc2      = pc_mean2(:,:,:,1);

% identify and ignore nans
nan_idx_cnt1   = ~isnan(d_behav1); 
nan_idx_cnt2   = ~isnan(d_behav2);
nan_idx_meg1   = ~any(isnan(squeeze(pc_mean1(1,2,:,:))),2);
nan_idx_meg2   = ~any(isnan(squeeze(pc_mean2(1,2,:,:))),2);
nan_idx1       = nan_idx_cnt1(:)&nan_idx_meg1(:);
nan_idx2       = nan_idx_cnt2(:)&nan_idx_meg2(:);

d_behav1    = permute(repmat(d_behav1(:),[1 400 400]),[2 3 1]);
d_behav2    = permute(repmat(d_behav2(:),[1 400 400]),[2 3 1]);
d_behav    = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);


[r_cnt,p_cnt] = tp_corr(d_fc(:,:,nan_idx1),d_behav(:,:,nan_idx1),3);
[r_cnt1,p_cnt1] = tp_corr(d_fc1(:,:,nan_idx1),d_behav1(:,:,nan_idx1),3);
[r_cnt2,p_cnt2] = tp_corr(d_fc2(:,:,nan_idx1),d_behav2(:,:,nan_idx2),3);


mask      = logical(tril(ones(size(d_fc1,1),size(d_fc1,1)),-1));
% r12(ifoi) = corr(r_cnt1(mask),r_cnt2(mask))

figure; set(gcf,'color','w');

subplot(3,2,1); imagesc(r_cnt,[-0.3 0.3]); axis square; title('Counting: Average')
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(3,2,2); imagesc(r_cnt.*(p_cnt<0.05),[-0.3 0.3]); axis square;
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

subplot(3,2,3); imagesc(r_cnt1,[-0.3 0.3]); axis square; title(sprintf('Counting: Block #1 (f = %d',ifoi))
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5)
subplot(3,2,4); imagesc(r_cnt1.*(p_cnt1<0.05),[-0.3 0.3]); axis square;
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

subplot(3,2,5); imagesc(r_cnt2,[-0.3 0.3]); axis square;  title(sprintf('Counting: Block #2 (f = %d)',ifoi))
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; set(gca,'FontSize',5);
subplot(3,2,6); imagesc(r_cnt2.*(p_cnt2<0.05),[-0.3 0.3]); axis square;
% set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
% set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; set(gca,'FontSize',5); tp_colorbar('Correlation');

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_behav_placebo_cortex_blocks_f%d_v%d.pdf',ifoi,v))
end

