%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 12;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

error('!')

%%

s = squeeze(nanmean(nanmean(squeeze(cleandat(:,:,:,2,2,6)-cleandat(:,:,:,1,2,6)))));
s1 = squeeze(nanmean(squeeze(cleandat(:,:,:,2,2,6)-cleandat(:,:,:,1,2,6))));

para.cond = 'cnt';
para.subj = SUBJLIST;

behav = pconn_getbehavior(para);

d = behav(:,2)-behav(:,1);

for i = 1 : 400
  i
    
    [r(i),p(i)] = corr(squeeze(s1(i,:)'),squeeze(d));
    
end

%%








