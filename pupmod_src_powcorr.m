%% pupmod_src_powcorr

clear all

% --------------------------------------------------------
% VERSION 1 - alpha0 = 0.05 
% --------------------------------------------------------
v         = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid      = 'cortex_lowres';
foi_range = 2.^(2:.25:6);
REG       = 0.05;
% --------------------------------------------------------
% VERSION 2 - alpha0 = 0.15
% --------------------------------------------------------
% v         = 2;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 0.15;
% --------------------------------------------------------
% VERSION 3 - alpha0 = 0.3 (version in main figures of the paper)
% --------------------------------------------------------
% v         = 3;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 0.3;
% --------------------------------------------------------
% VERSION 4 - alpha0 = 1.0
% --------------------------------------------------------
% v         = 4;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 1;
% --------------------------------------------------------
% VERSION 33 - AAL and alpha0 = 0.3
% --------------------------------------------------------
% v         = 33;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'aal_6mm';
% foi_range = 2.^(2.75:.25:4.25); % freq range based on sign. effects 
% REG       = 0.3;
% --------------------------------------------------------

if strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
end

outdir   = '/home/tpfeffer/pupmod/proc/conn/';

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % ------------
    % Create empty text file for parallelizatiopn
    % ------------
    if ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_src_powcorr_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
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
      
      % ------------
      % Load source model
      % ------------
      sa        = load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      % ------------

      for ifoi = 1:length(foi_range)
        
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
        [powcorr(:,:,iblock,ifoi), variance(:,iblock,ifoi)] = tp_data2orthopowcorr_wavelet(dat,filt,para);
        % ------------
        clear cs para filt

      end 
    end
    
    save(sprintf([outdir 'pupmod_src_variance_s%d_m%d_v%d.mat'],isubj,m,v),'variance');
    save(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_v%d.mat'],isubj,m,v),'powcorr');
    clear powcorr        
    
  end
end
error('!')


