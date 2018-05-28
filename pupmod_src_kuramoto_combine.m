
%% COMBINE RESULTS FROM TASK / REST
% Kuramoto order parameter computed in 
% pupmod_src_kuramoto.m
% pupmod_task_src_kuramoto.m

v = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

ord = pconn_randomization;
addpath ~/pconn/matlab
outdir = '~/pupmod/proc/conn/';

for isubj = SUBJLIST
  isubj
  for m = 1:3
    im = find(ord(isubj,:)==m);
    for iblock = 1 :2
      for ifoi = 1:2
      
        try    
          load(sprintf([outdir 'pupmod_src_kuramoto_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
          kura_std(isubj,m,iblock,ifoi,1) = std(R);
          kura_mean(isubj,m,iblock,ifoi,1) = mean(R);
        catch me
          kura_std(isubj,m,iblock,ifoi,1) = nan;
          kura_mean(isubj,m,iblock,ifoi,1) = nan;
        end

        clear R


        try    
          load(sprintf([outdir 'pupmod_task_src_kuramoto_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
          kura_std(isubj,m,iblock,ifoi,2) = std(R);
          kura_mean(isubj,m,iblock,ifoi,2) = mean(R);
        catch me
          kura_std(isubj,m,iblock,ifoi,2) = nan;
          kura_mean(isubj,m,iblock,ifoi,2) = nan;
        end

        clear R
      
      end
      
    end
  end
end

kura_std  = squeeze(nanmean(kura_std(SUBJLIST,:,:,:,:),3));
kura_mean = squeeze(nanmean(kura_mean(SUBJLIST,:,:,:,:),3));

save(sprintf(['~/pupmod/proc/conn/' 'pupmod_all_kuramoto_v%d.mat'],v),'kura_std','kura_mean');


      
      
  
