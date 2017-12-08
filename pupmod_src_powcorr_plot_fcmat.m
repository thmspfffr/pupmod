%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 1;
restoredefaultpath

addpath ~/pconn/matlab
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
% addpath /home/gnolte/neuconn/OLD/matlab/rest/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

grid  = sa_meg_template.grid_cortex3000;
g1    = sa_meg_template.grid_cortex3000;
g2    = sa_meg_template.cortex10K.vc;
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = sa_meg_template.grid_cortex3000;


SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';
  
ord = pconn_randomization;

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v))

s_fc = cleandat; clear cleandat

%% PLOT FC MATRICES
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap1 = cbrewer('seq', 'Blues', 100,'pchip');% colormap(autumn);
cmap2 = cbrewer('seq', 'Oranges', 100,'pchip');% colormap(autumn)
cmap = [cmap1(end:-1:1,:); ones(50,3); cmap2];

clim = [-0.01 0.01];

ifoi = 7;

fc_atx(:,:,1) = nanmean(s_fc(:,:,:,2,1,ifoi),3)-nanmean(s_fc(:,:,:,1,1,ifoi),3);
fc_atx(:,:,2) = nanmean(s_fc(:,:,:,2,2,ifoi),3)-nanmean(s_fc(:,:,:,1,2,ifoi),3);
fc_atx_dd = fc_atx(:,:,2)-fc_atx(:,:,1);

fc_dpz(:,:,1) = nanmean(s_fc(:,:,:,3,1,ifoi),3)-nanmean(s_fc(:,:,:,1,1,ifoi),3);
fc_dpz(:,:,2) = nanmean(s_fc(:,:,:,3,2,ifoi),3)-nanmean(s_fc(:,:,:,1,2,ifoi),3);
fc_dpz_dd = fc_dpz(:,:,2)-fc_dpz(:,:,1);

fc_atx_thresh(:,:,1) = ttest(atanh(s_fc(:,:,:,2,1,ifoi)),atanh(s_fc(:,:,:,1,1,ifoi)),'dim',3);
fc_atx_thresh(:,:,2) = ttest(atanh(s_fc(:,:,:,2,2,ifoi)),atanh(s_fc(:,:,:,1,2,ifoi)),'dim',3);

fc_dpz_thresh(:,:,1) = ttest(atanh(s_fc(:,:,:,3,1,ifoi)),atanh(s_fc(:,:,:,1,1,ifoi)),'dim',3);
fc_dpz_thresh(:,:,2) = ttest(atanh(s_fc(:,:,:,3,2,ifoi)),atanh(s_fc(:,:,:,1,2,ifoi)),'dim',3);

fc_atx_dd_thresh = ttest(s_fc(:,:,:,2,2,ifoi)-s_fc(:,:,:,1,2,ifoi),s_fc(:,:,:,2,1,ifoi)-s_fc(:,:,:,1,1,ifoi),'dim',3);
fc_dpz_dd_thresh = ttest(s_fc(:,:,:,3,2,ifoi)-s_fc(:,:,:,1,2,ifoi),s_fc(:,:,:,3,1,ifoi)-s_fc(:,:,:,1,1,ifoi),'dim',3);


% WITH THRESHOLD

figure;

subplot(2,2,1)
imagesc(triu(fc_atx(:,:,1).*fc_atx_thresh(:,:,1),1)+rot90(rot90(triu(fc_atx(:,:,2).*fc_atx_thresh(:,:,2),1))),clim);
axis square off; colormap(cmap)

subplot(2,2,2)
imagesc(triu(fc_dpz(:,:,1).*fc_dpz_thresh(:,:,1),1)+rot90(rot90(triu(fc_dpz(:,:,2).*fc_dpz_thresh(:,:,2),1))),clim);
axis square off; colormap(cmap)

subplot(2,2,3)
imagesc(triu(fc_atx_dd.*fc_atx_dd_thresh,1),clim);
axis square off; colormap(cmap)

subplot(2,2,4)
imagesc(triu(fc_dpz_dd.*fc_dpz_dd_thresh,1),clim);
axis square off; colormap(cmap)

print(gcf,'-dpdf',sprintf('~/pupmod_fc_withthresh_f%d_v%d.pdf',ifoi,v))


% NO THRESHOLD
figure;

subplot(2,2,1)
imagesc(triu(fc_atx(:,:,1),1)+rot90(rot90(triu(fc_atx(:,:,2),1))),clim);
axis square off; colormap(cmap)

subplot(2,2,2)
imagesc(triu(fc_dpz(:,:,1),1)+rot90(rot90(triu(fc_dpz(:,:,2),1))),clim);
axis square off; colormap(cmap)

subplot(2,2,3)
imagesc(triu(fc_atx_dd,1),clim);
axis square off; colormap(cmap)

subplot(2,2,4)
imagesc(triu(fc_dpz_dd,1),clim);
axis square off; colormap(cmap)

print(gcf,'-dpdf',sprintf('~/pupmod_fc_nothresh_f%d_v%d.pdf',ifoi,v))
  
%% PERMUTATION TEST

nperm = 1000;


dat_cnt1 = s_fc(:,:,:,[1 2],2,:);
dat_res1 = s_fc(:,:,:,[1 2],1,:);
dat_cnt2 = s_fc(:,:,:,[1 3],2,:); clear s_cnt
dat_res2 = s_fc(:,:,:,[1 3],1,:); clear s_res

taskvsrest(:,:,:,1,:) = dat_res1(:,:,:,1,:);
taskvsrest(:,:,:,2,:) = dat_cnt1(:,:,:,1,:);

for iperm = 1 : nperm
  iperm
    all_idx1 = randi(2,[size(SUBJLIST,2),1]);

    idx1 = all_idx1;
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
      permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);
      
      permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
      permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:);
      
      permdat_tvr(:,:,i,1,:) = taskvsrest(:,:,i,idx1(i),:);
      permdat_tvr(:,:,i,2,:) = taskvsrest(:,:,i,idx2(i),:);

      
    end
    
    [h,~,~,s] = ttest(permdat_res2(:,:,:,2,7),permdat_res2(:,:,:,1,7),'dim',3);
    max_cnt1_atx(iperm) = max(max(abs(triu(s.tstat))));
    
end






