%% pupmod_src_dfa

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
% foi_range       = [11 16];
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
% VERSION 9 - CORTEX LOW RES
% --------------------------------------------------------
v               = 9;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = [11 16];
para.segleng    = 9 ./ foi_range;
para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
para.epleng     = 60;
lpc             = 0;
timevariant     = 0;
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 1;
allpara.tau     = nan;
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
elseif strcmp(allpara.grid,'cortex_lowres')
  v_grid = 9;
end

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    for ifoi = 1:length(foi_range)
% %         
      if ~exist(sprintf([outdir 'pupmod_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
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
        
        if strcmp(allpara.grid,'cortex_lowres')
          sa.sa.grid_cortex_lowres = select_chans(sa.sa.grid_cortex3000,400);
        end
        
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
        dfa = tp_dfa_ortho_weight(dat,pars,sa);
        
        if size(dfa,1) < 100 && size(dfa,1) > 80
          pars = [];
          pars.grid = 'medium';
          dfa = tp_match_aal(pars,dfa);
        end
        
        dfa = dfa.exp;
%    
        save(sprintf([outdir 'pupmod_src_dfa_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'dfa');
        
      end
    end
  end
end
  
  
error('!')

%% 

clear dfa_all;

ord = pconn_randomization;
for isubj = SUBJLIST
  for m = 1:3
    im = find(ord(isubj,:)==m);
    for ifoi = 1:length(foi_range)

      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
      for iblock = 1:2
        
        load(sprintf([outdir 'pupmod_src_dfa_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));

        dfa_all(:,isubj,m,ifoi,iblock) = dfa(1:90);
        
      end
    end
  end
end

dfa_all = dfa_all(:,SUBJLIST,:,:,:);

%%

ttest(squeeze(mean(mean(dfa_all(:,:,1,1),5),1)),squeeze(mean(mean(dfa_all(:,:,2,1),5),1))
)
