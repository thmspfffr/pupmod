v= 1;
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
outdir = '~/pupmod/proc/sens/'
addpath ~/pconn/matlab/
ord   = pconn_randomization;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
     if ~exist(sprintf([outdir 'pupmod_sens_dfa_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_sens_dfa_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
  im = find(ord(isubj,:)==m);
  for iblock = 1 : 2
    
   
    
    load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))

    siginfo = nbt_Info;
    siginfo.converted_sample_frequency = 400;
          
          % compute bp-filtered signal
    tmp       = single(nbt_filter_fir(dat(:,~isnan(dat(1,:)))',8,12,400,2/8));
          ampenv    = abs(hilbert(tmp)); clear tmp mydata
          
          tmp_dfa   = tp_dfa(ampenv,[3 50],400,.5,15);
    
    dfa(iblock) = mean(tmp_dfa.exp);
    
    
  end
  save(sprintf('~/pupmod/proc/sens/pupmod_sens_dfa_s%d_m%d_v%d.mat',isubj,m,v),'dfa')
  end
end

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
     if ~exist(sprintf([outdir 'pupmod_cnt_sens_dfa_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_cnt_sens_dfa_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
  for iblock = 1 : 2
    
   
    
    load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))

    siginfo = nbt_Info;
    siginfo.converted_sample_frequency = 400;
          
          % compute bp-filtered signal
    tmp       = single(nbt_filter_fir(dat(:,~isnan(dat(1,:)))',8,12,400,2/8));
          ampenv    = abs(hilbert(tmp)); clear tmp mydata
          
          tmp_dfa   = tp_dfa(ampenv,[3 50],400,.5,15);
    
    dfa(iblock) = mean(tmp_dfa.exp);
    
    
  end
  save(sprintf('~/pupmod/proc/sens/pupmod_cnt_sens_dfa_s%d_m%d_v%d.mat',isubj,m,v),'dfa')
  end
end

 %%
 for isubj = SUBJLIST
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    load(sprintf('~/pupmod/proc/sens/pupmod_sens_dfa_s%d_m%d_v%d.mat',isubj,im,v))  
    dfa_all(isubj,m,1,:)=dfa; clear dfa


    load(sprintf('~/pupmod/proc/sens/pupmod_cnt_sens_dfa_s%d_m%d_v%d.mat',isubj,m,v))
    dfa_all(isubj,m,2,:)=dfa;  clear dfa

 
  end
 end