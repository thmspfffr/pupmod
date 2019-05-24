% pupmod_sens_peakfreq


clear

v = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
outdir = '~/pupmod/proc/sens/';

for isubj = SUBJLIST
  for m = 1:3
    
    if ~exist(sprintf([outdir 'pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    fprintf('Processing s%d m%d f%d ...\n', isubj,m)
    
    for iblock = 1:2
      
      fprintf('Loading MEG data ...\n');
      
      
      try
        load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        %           out.pf =
        %           save(sprintf([outdir 'pupmod_sens_peakfreq_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        continue
      end
      
      dat = data.trial';
      
      if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
        if (data.end_of_recording-600*data.fsample)<1
          data.start_of_recording = 1;
        else
          data.start_of_recording = data.end_of_recording-600*data.fsample;
        end
      elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
        if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
          data.end_of_recording = size(data.trial,2);
        else
          data.end_of_recording = data.start_of_recording+600*data.fsample;
        end
      elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      else
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      end
      
      dat = dat(data.start_of_recording:data.end_of_recording,:);
      dat(isnan(dat(:,1)),:)=[];
      
      for ichan = 1 : size(dat,2)
        %           ichan
        dx = abs(fft(dat(:,ichan))).^2; dx = dx(1:length(dx)/2+1);
        para.f = 0:400/size(dat,1):400/2;
        para.f_select = para.f >7 & para.f <13;
        para.method = 'gaussian';
        para.detrend = 0;
        para.detrendfreq = [1 100];
        [out.pf(ichan) out.fun(:,ichan)] = tp_peakfreq(dx,para);
        
      end
      
      dx = mean(abs(fft(dat)).^2,2);
      [out.pf_avg out.fun_avg] = tp_peakfreq(dx,para);
      
      save(sprintf('~/pupmod/proc/sens/pupmod_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v),'out');
      
      clear out
      
    end
  end
end

%%
addpath ~/pconn/matlab
clear pf_all pf_avg_all
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

v = 1
ord = pconn_randomization;
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      im = find(ord(isubj,:)==m);
      try
        load(sprintf('~/pupmod/proc/sens/pupmod_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
        
        if size(out.pf,2)==274
          pf_all(:,isubj,m,iblock) = out.pf;
        else
          pf_all(:,isubj,m,iblock) = pconn_sens_interp274_alena(isubj,out.pf,im);
        end
        pf_avg_all(isubj,m,iblock) = out.pf_avg;
      catch me
        
        pf_all(:,isubj,m,iblock) = nan(274,1);
        pf_avg_all(isubj,m,iblock)  = nan;
      end
    end
    
    
  end
end

pf_all = pf_all(:,SUBJLIST,:,:);
pf_avg_all= pf_avg_all(SUBJLIST,:,:);

%%

para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt = behav;

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);

%% correlate drug effects
ipharm = 2;
iblock = 1:2;
ifoi = 6;

tmp_fc = nanmean(nanmean(fc(:,:,:,:,2,ifoi,iblock),7),6);

d_beh = nanmean(behav_cnt(ipharm,:,iblock)-behav_cnt(1,:,iblock),3);
d_pf = nanmean(b(:,ipharm,iblock)-b(:,1,iblock),3);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7)-nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,2,ifoi,iblock),1),2),7),6));
% 
% d_beh = (nanmean(behav_cnt(ipharm,:,iblock)-behav_cnt(1,:,iblock),3))./nanmean(behav_cnt(1,:,iblock),3);
% d_pf = (nanmean(b(:,ipharm,iblock)-b(:,1,iblock),3))./nanmean(b(:,1,iblock),3);
% d_meg = (nanmean(reshape(tmp_fc(:,:,:,ipharm),[400*400 28]))-nanmean(reshape(tmp_fc(:,:,:,1),[400*400 28])))./nanmean(reshape(tmp_fc(:,:,:,1),[400*400 28]));


nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

figure; set(gcf,'color','w');
subplot(1,2,1);
scatter(d_beh,d_pf,50,'markerfacecolor','k','markeredgecolor','w'); 
xlabel('\DeltaSwitches'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_beh(nan_idx)',d_pf(nan_idx))
text(-38,-0.7,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)


subplot(1,2,2);
scatter(d_meg,d_pf,50,'markerfacecolor','k','markeredgecolor','w'); 
xlabel('\DeltaFC'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_meg(nan_idx),d_pf(nan_idx))
text(-0.018,-0.7,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)


%% correlate baseline
ipharm = 2;
iblock = 1:2;

ifoi = 6;

d_beh = nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1);
d_pf = nanmean(nanmean(b(:,ipharm,iblock),3),2);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7),4));

nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

figure; set(gcf,'color','w');
subplot(1,2,1);
scatter(d_beh,d_pf); 
xlabel('\DeltaSwitches'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_beh(nan_idx)',d_pf(nan_idx))
text(2,8.8,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)

subplot(1,2,2);
scatter(d_meg,d_pf); 
xlabel('\DeltaFC'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_meg(nan_idx),d_pf(nan_idx))
text(0.001,8.8,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)

%% topos

ipharm = 2;
iblock = 1:2;

ifoi = 6;

d_beh = nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1);
d_beh = repmat(d_beh,[274,1]);
d_pf = nanmean(nanmean(pf_all(:,:,ipharm,iblock),4),3);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7),4));
d_meg = repmat(d_meg',[274 1]);

nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

dd_pf = 100*(nanmean(nanmean(pf_all(:,:,ipharm,iblock),4),3)-nanmean(nanmean(pf_all(:,:,1,iblock),4),3))./nanmean(nanmean(pf_all(:,:,1,iblock),4),3);


% load /home/tpfeffer/pconn/proc/src/pconn_sa_s34_m1_b1_v1.mat
[r,p]=corr(d_beh',d_pf');

pars = [];
figure; set(gcf,'color','white');
pars.scale=[-5 5];
pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.cmap = cmap;
pars.resolution = 300;
showfield_colormap(mean(dd_pf,2),sa.locs_2D,pars);