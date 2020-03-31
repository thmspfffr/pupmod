%% pupmod_src_variance_task

clear

% --------------------------------------------------------
% VERSION 23 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v               = 23;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = 2.^[1:.25:7];
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
width           = 4;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 400;

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid, 'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'genemaps')
  v_grid = 13;
elseif strcmp(allpara.grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(allpara.grid,'cortex800')
  v_grid = 16;
end

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults()

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % %
    if ~exist(sprintf([outdir 'pupmod_src_power_variance_task_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_src_power_variance_task_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    fprintf('Processing s%d m%d f%d ...\n', isubj,m)
    
    for iblock = 1:2
      
      pars = [];
      pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      sa        = load(pars.sa);
      
      fprintf('Loading MEG data ...\n');
      
      try
        load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
          if strcmp(allpara.grid,'cortex_lowres')
            outp.var(:,:,iblock) = nan(400,length(foi_range));
          end
          continue
        else
          error('Data corrupt?')
        end
      end
      
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
      
      data.trial = data.trial(:,data.start_of_recording:data.end_of_recording);
      data.time = data.time(data.start_of_recording:data.end_of_recording);      
      
      for ifoi = 1:length(foi_range)
        
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)
        
        data.time(isnan(data.trial(1,:)))=[];
        data.trial(:,isnan(data.trial(1,:)))=[];
        
        cs = tp_wavelet_crossspec(data,foi_range(ifoi),width);
        
        para.iscs = 1;
        para.reg  = 0.05;
        
        sa.sa.filt      = pconn_beamformer(real(cs),sa.sa.L_coarse,para);
        fprintf('Computing pow corr  ...\n')
        
        f = foi_range(ifoi);
        
        KERNEL = tp_mkwavelet(f,0.5,data.fsample);
        
        n_win = size(KERNEL,1);
        n_shift = round(0.5*n_win);
        nseg = floor((size(data.trial,2)-n_win)/n_shift+1);
        
        clear datasf1
        
        for j=1:nseg
          dloc2=data.trial(:,(j-1)*n_shift+1:(j-1)*n_shift+n_win)';
          if any(isnan(dloc2(:,1)))
            warning('NaN detected')
            continue
          end
          dataf=dloc2'*KERNEL;
          datasf1(:,j)=dataf'*sa.sa.filt;
        end
        
        outp.var(:,ifoi,iblock) = var(abs(datasf1).^2,[],2);
        outp.pow(:,ifoi,iblock) = mean(abs(datasf1).^2,2);

      end
    end
    %
    save(sprintf([outdir 'pupmod_src_power_variance_task_s%d_m%d_v%d.mat'],isubj,m,v),'outp');
    
  end
end

error('!')

%% CLEAN NON PROCESSED FILES
outdir   = '/home/tpfeffer/pupmod/proc/conn/';

cnt = 0;
v = 6;
cnt_exist = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        ifoi
        if exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v)) && exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt_exist = cnt_exist + 1;
          
          continue
        elseif exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        elseif exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
          delete(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt = cnt + 1;
        elseif ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          warning('Nothing exists')
          cnt = cnt+1;
        end
      end
    end
  end
end
cnt


