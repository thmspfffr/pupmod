%% PLOT PUPIL AND BEHAVIOR

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

% v  = 23;
% cleandat = pupmod_loadpowcorr(v,SUBJLIST,0);
% fc = cat(4,fc(:,:,:,:,:,:,1),fc(:,:,:,:,:,:,2));
% fc = (fc-nanmean(fc,4))./nanstd(fc,[],4);
% fc = cat(7,fc(:,:,:,1:3,:,:),fc(:,:,:,4:6,:,:));
%
para.str_behav  = 'count';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt       = behav;

% tmp = reshape(behav_cnt,[3*2 28]);
% tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
% behav_cnt = reshape(tmp,[3 28 2]);

para.str_behav  = 'numb_switches';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn      = behav;
behav_bttn      = permute(behav_bttn,[2 1 3]);
%
behav_pooled = (behav_cnt+behav_bttn)./2;
behav_pooled(isnan(behav_bttn))=behav_cnt(isnan(behav_bttn));
behav_pooled(isnan(behav_cnt))=behav_bttn(isnan(behav_cnt));

para= [];
para.time = [4000 7000];
pup = pp_loadpupil(para);
% pup(:,:,:,1) = rest | :,2) = counting 

% zscore pupil within subjects (across blocks, sessions, separately for rest and task)
for isubj = 1 : size(pup,1)
  tmp(isubj,:,:,1) = reshape((reshape(pup(isubj,:,:,1),[6 1])-nanmean(reshape(pup(isubj,:,:,1),[6 1])))./nanstd(reshape(pup(isubj,:,:,1),[6 1])),[3 2]); 
  tmp(isubj,:,:,2) = reshape((reshape(pup(isubj,:,:,2),[6 1])-nanmean(reshape(pup(isubj,:,:,2),[6 1])))./nanstd(reshape(pup(isubj,:,:,2),[6 1])),[3 2]); 
end

clear pup; 
pup = tmp;

pup_pooled = nanmean(pup(:,:,:,1:2),4);
% pup_pooled(20,:,:)= [];

load ~/pupmod/paramsPhys_zscored.mat
% excluded subjects: 1,3,12  (excluded due to poor performance; see methods)
paramsPhys.pupilMean([1 3 12],:) = [];
paramsPhys.pupilVari([1 3 12],:) = [];
% exclude subjects with nans
paramsPhys.pupilMean(isnan(paramsPhys.pupilMean(:,1)),:)=[];
paramsPhys.pupilVari(isnan(paramsPhys.pupilVari(:,1)),:)=[];

pup_matching = paramsPhys.pupilMean;

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%%
figure; set(gcf,'color','w');
subplot(3,3,1:2); hold on

par_pooled = squeeze(nanmean(pup_pooled,3));

r=rand(28,1)-0.5;

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
end

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(0,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(2,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(4,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

para = [];
para.alpha = 0.05;
para.tail = 0;
para.nperm = 10000;

[~,p_pboatx_pooled] = tp_permtest(par_pooled(:,2),par_pooled(:,1),para);
[~,p_pbodpz_pooled] = tp_permtest(par_pooled(:,3),par_pooled(:,1),para);
[~,p_atxdpz_pooled] = tp_permtest(par_pooled(:,2),par_pooled(:,3),para);

axis([-2 25 -2 2]); tp_editplots

% ------------------------
% BEHAVIOR PERCEPTUAL CHOICE, POOLED (BUTTON/COUNTING)
% ------------------------
subplot(3,3,4:5); hold on

par_pooled = squeeze(nanmean(behav_pooled,3))';

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
end

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(0,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(2,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(4,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

axis([-02 25 0 150]); tp_editplots; 

% ------------------------
% VALUE BASED CHOICE: PUPIL
% --------------------------------
subplot(3,3,1:2); hold on

for isubj=1:size(pup_matching,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj)+8,pup_matching(isubj,2),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+10,pup_matching(isubj,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj)+8 r(isubj)+10],[pup_matching(isubj,2) pup_matching(isubj,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
end

line([8 10],[nanmean(pup_matching(:,2)) nanmean(pup_matching(:,1))],'color','k')
plot(8,nanmean(pup_matching(:,2)),'o','markerfacecolor','w','markeredgecolor','k','markersize',3)
plot(10,nanmean(pup_matching(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

[~,p_pboatx_valuechoice] = tp_permtest(pup_matching(:,2),pup_matching(:,1),para);

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_pupil_behavior.pdf'))

%% SEPARATELY FOR REST AND TASK BLOCKS


cols = cbrewer('seq', 'Reds', 28);

figure; set(gcf,'color','w');
subplot(3,3,1:2); hold on

col_rand = (rand-0.5)/2;

% ------------------------------------
% REST PUPIL FIRST
% ------------------------------------
pup_rest = nanmean(pup(:,:,:,1),3);
par_pooled = squeeze(nanmean(pup_rest,3));

r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
end

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(0,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(2,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(4,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

axis([-2 25 -2 2]); tp_editplots

para = [];
para.alpha = 0.05;
para.tail = 0;
para.nperm = 10000;

[~,p_pboatx_res] = tp_permtest(par_pooled(:,2),par_pooled(:,1),para);
[~,p_pbodpz_res] = tp_permtest(par_pooled(:,3),par_pooled(:,1),para);
[~,p_atxdpz_res] = tp_permtest(par_pooled(:,2),par_pooled(:,3),para);

% ------------------------------------
% COOUNTING PUPIL SECOND
% ------------------------------------

pup_cnt = nanmean(pup(:,:,:,2),3);
par_pooled = squeeze(nanmean(pup_cnt,3));

r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj)+6,par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+8,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+10,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj)+6 r(isubj)+8],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+8 r(isubj)+10],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
end

line([6 8],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([8 10],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(6,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(8,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(10,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

para = [];
para.alpha = 0.05;
para.tail = 0;
para.nperm = 10000;

[~,p_pboatx_cnt] = tp_permtest(par_pooled(:,2),par_pooled(:,1),para);
[~,p_pbodpz_cnt] = tp_permtest(par_pooled(:,3),par_pooled(:,1),para);
[~,p_atxdpz_cnt] = tp_permtest(par_pooled(:,2),par_pooled(:,3),para);

axis([-2 25 -2 2]); tp_editplots

subplot(3,3,4:5); hold on

par_pooled = squeeze(nanmean(permute(behav_cnt,[2 1 3]),3));
par_pooled = par_pooled(~any(isnan(par_pooled),2),:);

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
end

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(0,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(2,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(4,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

par_pooled = squeeze(nanmean(permute(behav,[2 1 3]),3))';
par_pooled = par_pooled(~any(isnan(par_pooled),2),:);

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj)+6,par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+8,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+10,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj)+6 r(isubj)+8],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+8 r(isubj)+10],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
end

line([6 8],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([8 10],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

plot(6,nanmean(par_pooled(:,2)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(8,nanmean(par_pooled(:,1)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')
plot(10,nanmean(par_pooled(:,3)),'o','markersize',3,'markeredgecolor','k','markerfacecolor','w')

axis([-02 25 0 150]); tp_editplots; %lsline


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_pupil_behavior_resttaskseparately.pdf'))


%%

