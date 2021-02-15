% pupmod_sens_peakfreq


clear
% ---------
% VERSION 1: gaussian fit
% ---------
v = 1
method = 'gaussian';
detrend = 0;
% ---------
% VERSION 2: center of mass
% ---------
% v = 2;
% method = 'com';
% detrend = 0;
% ---------
% VERSION 3: center of mass
% ---------
% v = 3;
% method = 'gaussian';
% detrend = 1;
% ---------

addpath ~/Documents/MATLAB/fieldtrip-20160919/; ft_defaults

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
outdir = '~/pupmod/proc/sens/';
ord    = pconn_randomization;
%%
for isubj = SUBJLIST
  
  isubj
  
  m = find(ord(isubj,:)==1);
  
  if ~exist(sprintf([outdir 'pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt'],isubj,m,v))
    system(['touch ' outdir sprintf('pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
  else
    continue
  end
  
  for iblock =1:2
    
    
    % ------------
    % Load sensor level data
    % ------------
    try
      load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
    catch me
      if ~exist(sprintf('~/pupmod/proc/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        if strcmp(grid,'cortex_lowres')
          powcorr(:,:,iblock,1:length(foi_range)) = nan(400,400,2,length(foi_range));
        elseif strcmp(grid,'aal_6mm')
          powcorr(:,:,iblock,1:length(foi_range)) = nan(90,90,2,length(foi_range));
        end
        continue
      else
        error('Data corrupt?')
      end
    end
    
    dat = dat(:,~isnan(dat(1,:)));
    
    [px,f]=pwelch(dat',hanning(4000),0.5,5:0.01:15,400,'power');
    out.pwelch(:,:,iblock) = px;
    out.pwelchfreq = f;
    
  end
  
  save(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_peakfreq_s%d_v%d.mat',isubj,v),'out')
  
  clear out dat
  
  for iblock =1:2
    
    % ------------
    % Load sensor level data
    % ------------
    try
      load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
    catch me
      if ~exist(sprintf('~/pupmod/proc/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        if strcmp(grid,'cortex_lowres')
          powcorr(:,:,iblock,1:length(foi_range)) = nan(400,400,1,length(foi_range));
        elseif strcmp(grid,'aal_6mm')
          powcorr(:,:,iblock,1:length(foi_range)) = nan(90,90,1,length(foi_range));
        end
        continue
      else
        error('Data corrupt?')
      end
    end
    
    
    dat = dat(:,~isnan(dat(1,:)));
    
    [px,f]=pwelch(dat',hanning(4000),0.5,5:0.01:15,400,'power');
    out.pwelch(:,:,iblock) = px;
    out.pwelchfreq = f;
    
  end
  
  save(sprintf('~/pupmod/proc/sens/pupmod_task_sens_peakfreq_s%d_v%d.mat',isubj,v),'out')
  clear out dat
  
end

%% PLOT


for isubj = SUBJLIST
      
  load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_peakfreq_s%d_v%d.mat',isubj,v),'out')
  
  for iblock = 1 : 2
    pxx(:,1,isubj,iblock) = smooth(nanmean(out.pwelch(:,:,iblock),2),20);
  end
    
  load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_peakfreq_s%d_v%d.mat',isubj,v),'out')

  for iblock = 1 : 2
    pxx(:,2,isubj,iblock) = smooth(nanmean(out.pwelch(:,:,iblock),2),20);
  end
  
end

pxx = nanmean(pxx(:,:,SUBJLIST,:),4);

for i = 1 : size(pxx,3)
  
  [~,imax(i,1)]=max(pxx(out.pwelchfreq>7 & out.pwelchfreq<13,1,i));
  [~,imax(i,2)]=max(pxx(out.pwelchfreq>7 & out.pwelchfreq<13,2,i));
  
end


%% PLOT AVERAGE POWER SPECTRUM (ACROSS SUBJECTS)
smo = 20;
figure; set(gcf,'color','w'); hold on

plot(f,smooth(nanmean(pow(:,:,1,1),2),smo))
plot(f,smooth(nanmean(pow(:,:,1,2),2),smo))

idx=out.pwelchfreq>7&out.pwelchfreq<14;
[~,i1]=max(smooth(mean(pow(idx,:,1,1),2),smo));
[~,i2]=max(smooth(mean(pow(idx,:,1,2),2),smo));
ff=out.pwelchfreq(idx);
pf_diff=ff(i1)-ff(i2);

nperm=10000;
all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);

for iperm = 1 : nperm
  
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    permpow(:,i,1) = pow(:,i,1,idx1(i));
    permpow(:,i,2) = pow(:,i,1,idx2(i));
  end
  
  idx=out.pwelchfreq>7&out.pwelchfreq<14;
  [~,i1]=max(smooth(mean(permpow(idx,:,1),2),smo));
  [~,i2]=max(smooth(mean(permpow(idx,:,2),2),smo));
  ff=out.pwelchfreq(idx);
  
  perm_pf(iperm,1)=ff(i1);
  perm_pf(iperm,2)=ff(i2);
  
end


perm_diff= perm_pf(:,1)-perm_pf(:,2);

p=sum((pf(1)-pf(2))< perm_diff)/nperm

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