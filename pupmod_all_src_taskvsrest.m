%% pupmod_task_src_powcorr

clear

v = 12;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';
  
ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
% 
        s_fc(:,:,isubj,m,1,ifoi,iblock) = powcorr;
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
        s_fc(:,:,isubj,m,2,ifoi,iblock) = powcorr;

      end
    end
  end
end

s_fc = nanmean(s_fc(:,:,SUBJLIST,:,:,:,:),7);
  
error('!')

%%

for ifoi = 1 : 13
  tmp = squeeze(s_fc(:,:,:,1,:,ifoi));
  
  [h,~,~,s] = ttest(tmp(:,:,:,2),tmp(:,:,:,1),'dim',3);
  
  n_p(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)>0))/4005;
  n_n(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)<0))/4005;
  
end

%% permutation test

nperm = 10000;

for iperm = 1 : nperm
  iperm
  clear permdat
  
  idx(:,1) = randi(2,[28 1]);
  idx(:,2) = 3-idx(:,1);
  
  for isubj = 1 : 28
    
    permdat(:,:,isubj,1,:) = squeeze(s_fc(:,:,isubj,1,idx(isubj,1),:));
    permdat(:,:,isubj,2,:) = squeeze(s_fc(:,:,isubj,1,idx(isubj,2),:));
 
  end
  
  for ifoi = 1 : 13
    
    tmp = squeeze(permdat(:,:,:,:,ifoi));
    [h,~,~,s] = ttest(tmp(:,:,:,2),tmp(:,:,:,1),'dim',3);
    
    n_p_perm(iperm,ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)>0))/4005;
    n_n_perm(iperm,ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)<0))/4005;
    
    
  end
  
  
end



  
  
  
  

  
  

