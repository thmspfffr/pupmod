%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 1;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

ord = pconn_randomization;
s_fc = single(zeros(90,90,34,3,2,13));

for ifoi = 1:13
  ifoi
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
        p1(:,:,iblock) = single(powcorr);
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
        p2(:,:,iblock) = single(powcorr);
        
        
      end
      
      s_fc(:,:,isubj,m,1,ifoi) = nanmean(p1,3);
      s_fc(:,:,isubj,m,2,ifoi) = nanmean(p2,3);
      
      clear p1 p2
      
      
    end
  end
end

s_fc = s_fc(:,:,SUBJLIST,:,:,:);

error('!')

%%

s = squeeze(nanmean(nanmean(squeeze(s_fc(:,:,:,3,2,6)-s_fc(:,:,:,1,2,6)))));

para.cond = 'cnt';
para.subj = SUBJLIST;

behav = pconn_getbehavior(para);

d = behav(:,2)-behav(:,1);

for i = 1 : 90
  i
  for j = 1 : 90
    
    [r(i,j),p(i,j)] = corr(squeeze(s(i,j,:)),squeeze(d));
    
  end
end

%%








