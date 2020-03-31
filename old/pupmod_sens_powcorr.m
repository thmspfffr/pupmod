%% pupmod_sens_powcorr

clear all

% --------------------------------------------------------
% VERSION 1 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v         = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
foi_range = 2.^(2:.25:6);
% --------------------------------------------------------

outdir   = '/home/tpfeffer/pupmod/proc/conn/';

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % ------------
    % Create empty text file for parallelizatiopn
    % ------------
    if ~exist(sprintf([outdir 'pupmod_sens_powcorr_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_sens_powcorr_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
% %     
    for iblock = 1:2
      
      clear dat sa
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

      for ifoi = 1:length(foi_range)
        
        % ------------
        % Compute orthogonalized power enevelope correlations 
        % ------------
        para            = [];
        para.fsample    = 400;
        para.freq       = foi_range(ifoi);
        powcorr(:,:,iblock,ifoi) = tp_data2orthopowcorr_wavelet_sens(dat,para);
        % ------------
        clear cs para filt

      end 
    end
    
    save(sprintf([outdir 'pupmod_sens_powcorr_s%d_m%d_v%d.mat'],isubj,m,v),'powcorr');
    clear powcorr        
    
  end
end
error('!')


