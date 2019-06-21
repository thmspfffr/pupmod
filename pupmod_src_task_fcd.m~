%% pupmod_task_src_powcorr

% Call pupmod_all_powcorr_periphereal.m next, in order to clean
% estimated FCs from artifacts. 

clear

% --------------------------------------------------------
% VERSION 1 - WEIGHTED AAL
% --------------------------------------------------------
v               = 1;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal_4mm';
foi_range       = unique(round(2.^[1:.5:7]));
para.segleng    = 9 ./ foi_range;
para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
para.epleng     = 60;
lpc             = 0;
timevariant     = 0;
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 1;
allpara.tau     = 0;
% --------------------------------------------------------
% VERSION 2 - SUM OVER ALL AAL VOXELS
% --------------------------------------------------------
% v               = 2;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'xcoarse';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = 0; % in seconds (0 = 1 sample)
% --------------------------------------------------------
% VERSION 10 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 10;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 12;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'cortex_lowres';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 13 - ELORETA
% --------------------------------------------------------
% v               = 13;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'eloreta';
% allpara.grid    = 'aal_4mm';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 1;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 14 - LCMV + moving average
% --------------------------------------------------------
% v               = 14;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 1;
% allpara.reg     = 0.05;
% allpara.weigh   = 1;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 15 - GENE EXPRESSION MAPS
% --------------------------------------------------------
% v               = 15;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'genemaps_aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% --------------------------------------------------------

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

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/


%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1:length(foi_range)
      
      if ~exist(sprintf([outdir 'pupmod_task_src_fcd_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_task_src_fcd_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      
      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
      for iblock = 1:2
        
        pars = [];
        pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        sa        = load(pars.sa);
        
        fprintf('Loading MEG data ...\n');
        
        try
          load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          if ~exist(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
            if strcmp(allpara.grid,'cortex_lowres')
              powcorr = nan(400,400);
            elseif strcmp(allpara.grid,'aal_4mm')
              powcorr = nan(91,91);
            elseif strcmp(allpara.grid,'aal_6mm')
              powcorr = nan(91,91);
            elseif strcmp(allpara.grid,'genemaps_aal')
              powcorr = nan(378,378);
            end
            
            save(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
            continue
          else
            error('Data corrupt?')
          end
        end
        
        dat = megdata2mydata(data); clear data
        
        pars = [];
        
        pars.fsample   = 400;
        
        if strcmp(para.wavelet,'bp_filt')
          pars.segleng   = round(para.segleng.*fsample);
          pars.segshift  = round(fsample*para.segleng/2);
        else
          pars.segleng   = round(para.segleng(ifoi).*fsample);
          pars.segshift  = round(fsample*para.segleng(ifoi)/2);
        end
        
        if ~any(size(foi_range)==1)
          pars.foi       = foi_range(ifoi,:);
        else
          pars.foi       = foi_range(ifoi);
        end
        
        pars.epleng    = size(dat,1);
        pars.epshift   = pars.epleng;
        pars.grid      = allpara.grid;
        pars.wavelet   = para.wavelet;
        pars.scnd_filt = para.scnd_filt;
        pars.filt      = allpara.filt;
        pars.tau       = allpara.tau;
        pars.reg       = allpara.reg;
        pars.bpfreq    = para.bpfreq(ifoi,:);
        pars.weigh     = allpara.weigh;
        
        % COMPUTE POWER CORRELATIONS
        % dat should be n x nchans
        pars.epleng   = para.epleng*pars.fsample;
        pars.epshift  = round(pars.epleng/16);
        
        if allpara.weigh == 0
          [powcorr] = tp_powcorr_ortho(dat,pars,sa);
        else
          [powcorr] = tp_powcorr_ortho_weight(dat,pars,sa);
        end
        
        
        save(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        
      end
    end
  end
end


error('!')

%%
outdir   = '/home/tpfeffer/pupmod/proc/conn/';
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 34];
cnt = 0;
v = 13;
cnt_exist = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        %         ifoi
        if exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v)) && exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt_exist = cnt_exist + 1;
          
          continue
        elseif exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
          
        elseif exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
          %           delete(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt = cnt + 1;
          fprintf('S%dm%df%db%d\n',isubj,m,ifoi,iblock)
          
        elseif ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && exist(sprintf([outdir 'pupmod_task_src_fpowcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          warning('Nothing exists')
          cnt = cnt+1;
        end
      end
      
    end
  end
end
cnt





%% PLOT FILTERED TIME COURSE

flp = 8;           % lowpass frequency of filter
fhi = 12;

para.ord = 4;
delt = 1/400;            % sampling interval
k=4;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);

dat = data.trial{1}(1cl:20,:);

dat1 = hilbert(filtfilt(bfilt,afilt,dat'))

for i = 1 : 20
  figure; set(gcf,'color','w'); hold on
  plot(real(dat1(50000:50800,i)),'color',[0.7 0.7 0.7],'linewidth',2); 
  plot(abs(dat1(50000:50800,i)),'r','linewidth',3); 
  tp_editplots; 
  set(gca,'xtick',[0 400 800],'xticklabels',[0 1 2]); xlabel('Time [s]')
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_timeseries_i%d.pdf',i))
end