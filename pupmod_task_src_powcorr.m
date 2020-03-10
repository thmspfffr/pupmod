%% pupmod_task_src_powcorr

clear all

% --------------------------------------------------------
% VERSION 1 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v         = 1;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG         = 0.05;
% --------------------------------------------------------
% VERSION 1 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v         = 11;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 0.1;
% --------------------------------------------------------
% VERSION 111 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v         = 111;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid      = 'cortex_lowres';
foi_range = 2.^(2:.25:6);
REG       = 1;
% --------------------------------------------------------
% VERSION 2 - VOXEL LEVEL, 400 samples cortex, but not filtered
% --------------------------------------------------------
% v         = 2;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'aal_6mm';
% foi_range = 2.^(2:.25:6);
% REG       = 0.05;
% --------------------------------------------------------

if strcmp(grid,'xcoarse')
  v_grid = 2;
elseif strcmp(grid,'cortex')
  v_grid = 3;
elseif strcmp(grid,'aal')
  v_grid = 4;
elseif strcmp(grid,'medium')
  v_grid = 5;
elseif strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
elseif strcmp(grid,'genemaps')
  v_grid = 13;
elseif strcmp(grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(grid,'cortex800')
  v_grid = 16;
end

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
% %     
    if ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_task_src_powcorr_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
%     
    for iblock = 1:2
      
      clear dat sa
      % ------------
      % Load sensor level data
      % ------------
      try
        load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('~/pupmod/proc/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))             
          if strcmp(grid,'cortex_lowres')
            powcorr(:,:,iblock,1:length(foi_range)) = nan(400,400,1,length(foi_range));
          elseif strcmp(grid,'aal_6mm')
            powcorr(:,:,iblock,1:length(foi_range)) = nan(90,90,1,length(foi_range));
          end
          continue
        else
          error('Data corrupt?')
        end
      end
      % ------------
      % Load source model
      % ------------
      sa        = load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      % ------------

      for ifoi = 1 : length(foi_range)
        
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)
        % ------------
        % Compute cross spectrum (using wavelets)
        % ------------
        para          = [];
        para.freq     = foi_range(ifoi);
        para.fsample  = 400;
        cs            = tp_compute_csd_wavelets(dat,para);
        % ------------
        % Compute spatial filter (LCMV)
        % ------------
        para          = [];
        para.reg      = REG;
        filt          = tp_beamformer(real(cs),eval(sprintf('sa.sa.L_%s',grid)),para);   
        % ------------
        % Compute orthogonalized power enevelope correlations 
        % ------------
        para            = [];
        para.fsample    = 400;
        para.freq       = foi_range(ifoi);
        powcorr(:,:,iblock,ifoi) = tp_data2orthopowcorr_wavelet(dat,filt,para);
        % ------------
        clear cs para pow filt 

      end 
    end
    
    save(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_v%d.mat'],isubj,m,v),'powcorr');
    clear powcorr        
    
  end
end
error('!')