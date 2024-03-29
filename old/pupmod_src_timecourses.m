%% pupmod_src_powcorr

clear 

% --------------------------------------------------------
% VERSION 1 - WEIGHTED AAL
% --------------------------------------------------------
% v               = 1;
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
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 1;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 2 - SUM OVER ALL AAL VOXELS
% --------------------------------------------------------
% v               = 2;
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
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 2;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 2 - M758
% --------------------------------------------------------
% v               = 3;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'm758_4mm';
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
% VERSION 5 - XCOARSE
% --------------------------------------------------------
% v               = 5;
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
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 6 - XCOARSE
% --------------------------------------------------------
% v               = 6;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'cortex';
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
% allpara.tau     = 10;
% --------------------------------------------------------
% VERSION 10 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 10;
% v_postproc      = 6;
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
% allpara.reg = 0.05;
% allpara.weigh     = 0;
% --------------------------------------------------------


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
% addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 400;

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
end

% t = license('test','signal_toolbox');
% if t
% %   continue
% else
%   error('Toolbox not available');
% end
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1 : length(foi_range)
% %         
      if ~exist(sprintf([outdir 'pupmod_src_timecourses_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_timecourses_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
      for iblock = 1:2
        
        fprintf('Loading MEG data ...\n');
        
        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
        [dat] = megdata2mydata(data); clear data
        
        pars      = [];
        pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        sa        = load(pars.sa);
        
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
        
        pars.reg       = allpara.reg;
        pars.bpfreq    = para.bpfreq(ifoi,:);
        
        epleng    = size(dat,1);
        epshift   = pars.epleng;
        segleng   = para.segleng;
        segshift  = para.segshift;
        f         = para.foi;
        scale     = sqrt(mean(mean(data.^2)));
        data      = data/scale;
        [n nchan] = size(data);


        if any(size(f)~=1)
          para.freqoi = [f(1) f(2)];
          flp = f(1);           % lowpass frequency of filter
          fhi = f(2);
        else
          para.freqoi = [f-1 f+1];
          flp = f-1;           % lowpass frequency of filter
          fhi = f+1;
        end


        delt = 1/fsample; k=4;                  
        fnq=1/(2*delt); Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
        [bfilt,afilt]=butter(k,Wn);

        dloc2    = hilbert(zscore(filtfilt(bfilt,afilt,dat))); 

       save(sprintf([outdir 'pupmod_src_timecourses_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        
      end
    end
  end
end
  
  
error('!')


  
%%

  
