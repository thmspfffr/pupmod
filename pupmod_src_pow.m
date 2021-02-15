%% pupmopd_src_pow
% compute power in AAL parcellation using beamforming

clear
restoredefaultpath

% -------------------------
% VERSION 33 - AAL and alpha0 = 0.3
% --------------------------------------------------------
v         = 33;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid      = 'aal_6mm';
foi_range = 2.^(2:(1/4):6); % freq range based on sign. effects
REG       = 0.3;
% --------------------------------------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pupmod/matlab/

ft_defaults
addpath ~/pconn/matlab/
outdir = '~/pupmod/proc/src/';

ord    = pconn_randomization;

if strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
end

%%
% -------------------------
for isubj = SUBJLIST
  
  for m = 1 : 3
    
    for iblock = 1:2
      %
      fn = sprintf('pupmod_src_pow_s%d_m%d_b%d_v%d',isubj,m,iblock,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
      %
      fprintf('Processing subj%d block%d ...\n',isubj,iblock);
      
      try
        % load cleaned meg data
        load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
      catch me
        continue
      end
      
      load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      lf = sa.L_aal_6mm;
      
      for ifreq=1:length(foi_range)
        ifreq
        
        % -------------------------------
        % compute CSD
        % -------------------------------
        para          = [];
        para.freq     = foi_range(ifreq);
        para.fsample  = 400;
        para.overlap  = 0.5;
        csd           = tp_compute_csd_wavelets(dat,para);
        % -------------------------------
        % beamforming
        % -------------------------------
        para          = [];
        para.reg      = 0.05;
        [~,outp.pow(:,ifreq)] = tp_beamformer(real(csd),lf,para);
        % --------------
        
        % assign locations to AAL regions
        for ilabel = 1 : max(sa.aal_label)
          idx = sa.aal_label==ilabel;
          outp.pow_aal(ilabel,ifreq) = mean(outp.pow(idx,ifreq),1);
        end
        
      end
      
      save([outdir fn '.mat'],'outp')
      tp_parallel(fn,outdir,0)
      
      clear outp
    end
  end
end


% DO SAME FOR TASK DATA
for isubj = SUBJLIST
  
  for m = 1 : 3
    
    for iblock = 1:2
      %
      fn = sprintf('pupmod_task_src_pow_s%d_m%d_b%d_v%d',isubj,m,iblock,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
      %
      fprintf('Processing subj%d block%d ...\n',isubj,iblock);
      
      try
        % load cleaned meg data
        load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
      catch me
        continue
      end
      
      load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      lf = sa.L_aal_6mm;
      
      for ifreq=1:length(foi_range)
        ifreq
        
        % -------------------------------
        % compute CSD
        % -------------------------------
        para          = [];
        para.freq     = foi_range(ifreq);
        para.fsample  = 400;
        para.overlap  = 0.5;
        csd           = tp_compute_csd_wavelets(dat,para);
        % -------------------------------
        % beamforming
        % -------------------------------
        para          = [];
        para.reg      = 0.05;
        [~,outp.pow(:,ifreq)] = tp_beamformer(real(csd),lf,para);
        % --------------
        
        % assign locations to AAL regions
        for ilabel = 1 : max(sa.aal_label)
          idx = sa.aal_label==ilabel;
          outp.pow_aal(ilabel,ifreq) = mean(outp.pow(idx,ifreq),1);
        end
        
      end
      
      save([outdir fn '.mat'],'outp')
      tp_parallel(fn,outdir,0)
      
      clear outp
    end
  end
end


exit

%%
clear all_pow

v = 33

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1:2
      
      try
      load(sprintf([outdir 'pupmod_src_pow_s%d_m%d_b%d_v%d'],isubj,m,iblock,v));
      all_pow(:,isubj,m,iblock,:,1) = outp.pow_aal;
      catch me
        all_pow(:,isubj,m,iblock,:,1) = nan(91,17);
      end
      
      try
      load(sprintf([outdir 'pupmod_task_src_pow_s%d_m%d_b%d_v%d'],isubj,m,iblock,v));
      all_pow(:,isubj,m,iblock,:,2) = outp.pow_aal;
      catch me
        all_pow(:,isubj,m,iblock,:,2) = nan(91,17);
      end
    end
  end
end

all_pow = squeeze(nanmean(all_pow(:,SUBJLIST,:,:,:,:),4));


k = 1 : 90;

% exclude subcortical regions
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;
[~,trans,lab]=tp_match_aal(para,rand(90,90));

all_pow = all_pow(trans(:,2),:,:,:,:,:);
all_pow = all_pow(include_bcn,:,:,:,:,:);

[h,~,~,s] = ttest(nanmean(all_pow(:,:,1,6:9,2),4),nanmean(all_pow(:,:,1,6:9,1),4),'dim',2);

task_idx = h;

save(sprintf('~/pupmod/proc/src/pupmod_src_pow_taskmodulationindex_v%d.mat',v),'task_idx')


