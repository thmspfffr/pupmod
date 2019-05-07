%% CORRELATE PLACEBO BEHAVIOR WITH PLACEBO POWCORR


%% PLOT SPECTRUM OF CORRELATIONS WITH BEHAVIOR
mask      = logical(tril(ones(400,400),-1));

SUBJ = 1:28; %SUBJ(13)=[];
iblock = 1:2;
ipharm = 1:3;

for ifoi = 1:13
  ifoi
    
  d_behav    = nanmean(nanmean(behav_cnt(ipharm,SUBJ,iblock),3),1);
  d_behav = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);
  
  d_fc       = squeeze(nanmean(nanmean(fc(:,:,SUBJ,ipharm,1,ifoi,iblock),7),4));
  
  nan_idx_meg = ~isnan(squeeze(d_fc(1,2,:)));
  nan_idx_beh = ~isnan(squeeze(d_behav(1,2,:)));
  nan_idx     = nan_idx_meg&nan_idx_beh;

  % compute correlation
  [r,p]=tp_corr(d_fc(:,:,nan_idx),d_behav(:,:,nan_idx),3);
  pos_corr(ifoi) = 100* (sum(p(mask)<0.05 & r(mask)>0) / sum(mask(:)));
  neg_corr(ifoi) = 100* (sum(p(mask)<0.05 & r(mask)<0) / sum(mask(:)));
  
  pos_corr_vox(:,ifoi) = 100*(sum(p<0.05 & r>0)/sum(mask(:)));
  neg_corr_vox(:,ifoi) = 100*(sum(p<0.05 & r<0)/sum(mask(:)));
  
 pp(:,:,ifoi) = p;
 rr(:,:,ifoi) = r;

end
%%
figure; set(gcf,'color','w');

subplot(2,1,1); hold on

plot(pos_corr,'r');
plot(neg_corr,'b');

tp_editplots; 
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
xlabel('Carrier frequency [Hz]'); ylabel('Fraction of sign. correlations [%]')
%%
ifoi = 6;

para =[]
para.cmap = cmap;
% para.cmap = para.cmap(end:-1:1,:);
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.25 0.25];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_placebo_avg.png');
tp_plot_surface(nanmean(rr(:,:,ifoi))',para)


para =[]
para.cmap = plasma;
% para.cmap = para.cmap(end:-1:1,:);
para.grid = grid;
para.dd = 0.75;
para.clim = [0 0.1];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_placebo_avg.png');
tp_plot_surface(pos_corr_vox(:,ifoi)',para)


%%

