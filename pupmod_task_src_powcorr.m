%% pupmod_task_src_powcorr

clear all

% --------------------------------------------------------
% VERSION 1 - mean
% --------------------------------------------------------
% v               = 1;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.smo        = foi_range./4;
% para.segleng    = 1 ./ para.smo;
% para.epleng     = 5;
% lpc             = 0;
% timevariant     = 0;
% --------------------------------------------------------
% VERSION 2 - time variant
% --------------------------------------------------------
% v               = 2;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.smo        = foi_range./4;
% para.segleng    = 1 ./ para.smo;
% para.epleng     = 5;
% lpc             = 1;
% timevariant     = 0;
% para.wavelet    = '??'
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 3 - time variant
% --------------------------------------------------------
% v               = 3;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.smo        = foi_range./4;
% para.segleng    = 1 ./ para.smo;
% para.epleng     = 30;
% lpc             = 0;
% timevariant     = 1;
% para.wavelet    = '??'
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 4 - new wavelets
% --------------------------------------------------------
% v               = 4;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 5;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'ft';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 5 - new wavelets & lambda = 5% (instead of 1%)
% --------------------------------------------------------
% v               = 5;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 5;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'ft';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 6 - new wavelets & lambda = 5% (instead of 1%)
% --------------------------------------------------------
% v               = 6;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 5;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'ft';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 7 - new wavelets & lambda = 5% (instead of 1%)
% --------------------------------------------------------
% v               = 7;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 5;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'ft';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 8 - ft_wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 8;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 1;
% para.wavelet    = 'ft';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 9 - bandpass filter
% --------------------------------------------------------
% v               = 9;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 1;
% --------------------------------------------------------
% VERSION 10 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
v               = 10;
v_postproc      = 6;
v_grid          = 4; % 4 = aal
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal';
foi_range       = unique(round(2.^[1:.5:7]));
para.segleng    = 9 ./ foi_range;
para.epleng     = 60;
lpc             = 0;
timevariant     = 0;
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 11 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 11;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = [2 4; 4 8; 8 12; 12 24; 24 48; 52 100; 52 150];
% para.segleng    = 1;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 12 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 12;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'eloreta';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 13 - new wavelets & lambda = 5% (instead of 1%) - timevariant
% --------------------------------------------------------
% v               = 13;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'eloreta';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 1;
% --------------------------------------------------------
% VERSION 14 - same as v10, but with coherence
% --------------------------------------------------------
% v               = 14;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 15 - same as v14, just with zero-padding
% --------------------------------------------------------
% v               = 15;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% --------------------------------------------------------
% VERSION 16 - same as v14, just with zero-padding
% --------------------------------------------------------
% v               = 17;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 1;
% allpara.tau     = inf;
% --------------------------------------------------------
% VERSION 18 - same as v14, just with zero-padding
% --------------------------------------------------------
% v               = 18;
% v_postproc      = 6;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = [8 12];
% para.segleng    = 9 ./ foi_range;
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 1;
% allpara.tau     = inf;
% --------------------------------------------------------



addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 400;

t = license('test','signal_toolbox');
% if t
% %   continue
% else
%   error('Toolbox not available');
% end
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1 : length(foi_range)
%       
      if ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      
      for iblock = 1:2
        
        disp(sprintf('Loading MEG data ...'));
        
%         if isubj < 24
%           load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2));
%         else
%           load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
%         end
%          data_low.trial{1} = data_low.trial{1} + data_hi.trial{1}; clear data_hi
%         
%         [dat,epleng] = megdata2mydata(data_low); clear data_low
        
        try 
          load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          if ~exist(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
            powcorr = nan(90,90);
            save(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
            continue
          else
            error('Data corrupt?')
          end
        end
        
        [dat,epleng] = megdata2mydata(data); clear data
        
        pars = [];
        pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
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
        
        pars.epleng    = size(dat,1);
        pars.epshift   = pars.epleng;
        pars.aal       = 0;
        pars.wavelet   = para.wavelet;
        pars.scnd_filt = para.scnd_filt;
        pars.filt      = allpara.filt;
        pars.tau       = allpara.tau;

        % COMPUTE POWER CORRELATIONS
        % dat should be n x nchans
        if ~timevariant
          if ~lpc
            [powcorr] = tp_powcorr_ortho(dat,pars,sa);
          else
            powcorr = tp_data2lpc_jackknife(dat,pars,filt,filt);
          end
        else
          
          pars.epleng   = para.epleng*pars.fsample;
          pars.epshift  = round(pars.epleng/16);
          
          if ~lpc
            [powcorr] = tp_powcorr_ortho(dat,pars,sa);
          else
            powcorr = tp_lpc(dat,pars,filt,filt);
          end
        end

        pars = [];
        pars.grid = 'medium';
        powcorr = tp_match_aal(pars,powcorr);
%         coh = tp_match_aal(pars,coh);

       save(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        
      end
    end
  end
end
  
  
error('!')
  
%%
outdir   = '/home/tpfeffer/pupmod/proc/conn/';

cnt = 0;
v = 10;
cnt_exist = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        ifoi
        if exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v)) && exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt_exist = cnt_exist + 1;

          continue
        elseif exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);

        elseif exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
          delete(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt = cnt + 1;
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
  
%%
clear s s1 s2 fc_mean
v = 10;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

addpath ~/pconn/matlab/
  
ord = pconn_randomization;

for ifoi = 6:6
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,12,v));

        for i = 1 : size(powcorr,3)
          tmp(i)=nanmean(nanmean(powcorr(:,:,i)));
        end
        
        idx1 = find(tmp>median(tmp));
        idx2 = find(tmp<median(tmp));
        
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,5,v));
        
        s1(:,:,isubj,m,1,iblock) = mean(powcorr(:,:,idx1),3);
        s1(:,:,isubj,m,2,iblock) = mean(powcorr(:,:,idx2),3);

      end
    end
  end
end
%%

s = squeeze(nanmean(nanmean(nanmean(s1(:,:,SUBJLIST,:,:,:),6),3),4));
figure;
subplot(1,3,1)
imagesc(s(:,:,1),[-max(s(:)) max(s(:))]); axis square
subplot(1,3,2)
imagesc(s(:,:,2),[-max(s(:)) max(s(:))]); axis square
subplot(1,3,3)
imagesc(s(:,:,1)-s(:,:,2),[-0.01 0.01]); axis square
colormap(cmap)

% fc_mean = nanmean(fc_mean(SUBJLIST,:,1:13,:),4);

%% COMPUTE A PERMUTATION TEST FOR THIS

nfoi = 13;
clim1 = [-0.05 0.05];
nperm = 1000;

test_cond = 3;

% 
% for ifoi = 1 : nfoi
% 
%     t1=ttest(ss(:,:,:,1,ifoi),ss(:,:,:,2,ifoi),'dim',3);
%     t2=ttest(ss(:,:,:,1,ifoi),ss(:,:,:,3,ifoi),'dim',3);
% 
%     a1=t1.*(s2-s1); a1 = a1(:); p(ifoi,1) = sum(a1>0); n(ifoi,1) = sum(a1<0); 
%     a2=t2.*(s3-s1); a2 = a2(:); p(ifoi,2) = sum(a2>0); n(ifoi,2) = sum(a2<0); 
% 
% end



% figure; set(gcf,'color','w'); hold on

for iperm = 1 : nperm
    
  iperm
  
  idx1 = randi(2,[18,1]);
  idx2 = 3-idx1;
                 
  for i = 1 : length(idx1)
    
    ss(:,:,i,1,:) = s(:,:,i,idx1(i),:);
    ss(:,:,i,2,:) = s(:,:,i,idx2(i),:);
    
  end

  for ifoi = 1 : nfoi

    s1 = atanh(squeeze(nanmean(ss(:,:,:,1,ifoi),3)));
    s2 = atanh(squeeze(nanmean(ss(:,:,:,2,ifoi),3)));
%     s3 = atanh(squeeze(nanmean(ss(:,:,:,3,ifoi),3)));

    t1=ttest(ss(:,:,:,1,ifoi),ss(:,:,:,2,ifoi),'dim',3);
%     t2=ttest(ss(:,:,:,1,ifoi),ss(:,:,:,3,ifoi),'dim',3);

    a1=t1.*(s2-s1); a1 = a1(:); p_perm(ifoi,1,iperm) = sum(a1>0); n_perm(ifoi,1,iperm) = sum(a1<0); 
%     a2=t2.*(s3-s1); a2 = a2(:); p_perm(ifoi,2,iperm) = sum(a2>0); n_perm(ifoi,2,iperm) = sum(a2<0); 

%     clim1 = [min([min(s1(s1~=0)) min(s2(s2~=0))]) max([max(s1(s1~=0)) max(s2(s2~=0))])];
% 
%     subplot(nfoi,5,(ifoi-1)*5+1)
%     imagesc(s1,clim1); colormap(jet); colormap
%     subplot(nfoi,5,(ifoi-1)*5+2)
%     imagesc(s2,clim1); colormap(jet); colormap
%     subplot(nfoi,5,(ifoi-1)*5+3)
%     imagesc(s3,clim1); colormap(jet); colormap
%     subplot(nfoi,5,(ifoi-1)*5+4)
%     imagesc((s2-s1).*t1,[-0.01 0.01]); colormap(jet); colormap
%     subplot(nfoi,5,(ifoi-1)*5+5)
%     imagesc((s3-s1).*t2,[-0.01 0.01]); colormap(jet); colormap
  end
  clear ss idx1 idx2

end

% get p-values

for ifoi = 1 : 13
  
  pp(ifoi)=1-sum(p(ifoi,1)>p_perm(ifoi,1,:))./1000
  pn(ifoi)=1-sum(n(ifoi,1)>n_perm(ifoi,1,:))./1000

end



% figure; set(gcf,'color','white'); hold on;
% 
% plot(foi_range,p);
% plot(foi_range,n);
% 




%%
% LOAD SIGNIFICANT CLUSTER FROM DFA PROJECT 
% load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_dfa_c1_f%d_v%d.mat',3,1))
% load aalmask_grid_cortex3000
% 
% for i = 1 : 89
%   if sum(aalgrid.mask==i)==0
%     continue
%   end
%   
%   aalgrid.names(i) = unique(aalgrid.labels((aalgrid.mask==i)));
%   
%   a(i)=sum(stats.mask(aalgrid.mask==i))/sum(aalgrid.mask==i)>0.5;
% %   cnt = cnt + 1;
% end



a = logical(a);
s2 = s(a,:,:,2);
s1 = s(a,:,:,1);

nanmean(nanmean(squeeze(nanmean(s2,2))))
nanmean(nanmean(squeeze(nanmean(s1,2))))

% [t,p]=ttest(s2,s1,'dim',3);

% imagesc(t); 




%%
figure; set(gcf,'color','white');



clim = [0 0.25];

s_plt = nanmean(s.*1.77,6);
s_plt = s_plt(:,:,SUBJLIST,:,:);

ss = squeeze(nanmean(nanmean(nanmean(nanmean(s_plt,4),3),2),1));

%%











































%% PLOT PEAK FREQ OF VAR
a = squeeze(nanmean(s_plt,3));
a = squeeze(nanmean(a,3));
a = squeeze(nanmean(a,2));

for i = 1 : 91
%   for j = 1 : 91
    [~,t(i)]=max(a(i,:));
    idx(i) = foi_range(t(i));
%   end
end

%% PLOT ON AAL COARSE

load sa_meg_template;
mri   = sa_meg_template.mri;

addpath ~/pconn/matlab/
load aalmask_grid_cortex3000
grid = sa_meg_template.grid_cortex3000;

dat = zeros(3000,1);

for i = 1 : 91
  
  aalidx = find(aalgrid.mask==i);
  dat(aalidx) = idx(i);
  
end


g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

vc = sa_meg_template.vc;

viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

dd = .1;
%%
par_interp = spatfiltergauss(dat,g1,dd,g2);
   
figure; set(gcf,'color','white'); hold on;

para = [] ;
% r = max(abs(min(d)),abs(max(d)));
%   para.colorlimits = [min(d(d>eps)) max(d)];
para.colorlimits = [0 38];

for iplot = 1 : 6
  
  subplot(3,2,iplot)
  
  para.myviewdir = viewdir(iplot,:);
  a = sa_meg_template.cortex10K;
  
  if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
  elseif iplot == 6
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
  end
  
  pconn_showsurface(a,para,par_interp)
  
  colormap(parula)
  
  camlight headlight
  
end

% saveas(gcf,sprintf([plotdir 'pconn_src_contrast_f%d_v%d.fig'],ifoi,v),'fig')

%%













































%%

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)

ifoi = 6;


clim = [0 0.25];

s_plt = squeeze(nanmean(nanmean(nanmean(s(:,:,SUBJLIST,:,ifoi,:),6),5),3));

dat = zeros(3000,3);

sss = squeeze(nanmean(s_plt,2));

for i = 1 : 91
  
  aalidx = find(aalgrid.mask==i);
  dat(aalidx,:) = repmat(sss(i,:),[length(aalidx) 1]);
  
end

dat = dat(:,2)-dat(:,1);

g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

vc = sa_meg_template.vc;

viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

dd = .1;

par_interp = spatfiltergauss(dat,g1,dd,g2);
   
figure; set(gcf,'color','white'); hold on;

para = [] ;
% r = max(abs(min(d)),abs(max(d)));
%   para.colorlimits = [min(d(d>eps)) max(d)];
para.colorlimits = [-max([abs(min(dat)) abs(max(dat))]) max([abs(min(dat)) abs(max(dat))])];

for iplot = 1 : 6
  
  subplot(3,2,iplot)
  
  para.myviewdir = viewdir(iplot,:);
  a = sa_meg_template.cortex10K;
  
  if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
  elseif iplot == 6
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
  end
  
  pconn_showsurface(a,para,par_interp)
  
  colormap(cmap)
  
  camlight headlight
  
end



%%

for iblock = 1 : 2
  
  load(sprintf('/home/tpfeffer/pupmod/proc/pupmod_src_powcorr_s4_m1_b%d_f10_v1.mat',iblock))


  for i = 1  : 91
    for j = 1 : 91

      [r(i,j,iblock),p(i,j,iblock)]=corr(squeeze(par.powcorr(i,j,:)),par.pup');

    end
  end

end




  
  
  
  %% GET TIME COURSES IN ORDER TO COMPUTE LOCAL BIFURCATION PARAMETER
  % SOURCE SETTINGS
  para.grid = 'coarse';
  para.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
  para.filt = 'lcmv';
  para.fs   = 400;
  para.foi  = [8 12];
  para.cs   = cs;
  para.aal  = 0;
  
  [dat A1]  = pupmod_src_timecourse(mydata,para);
  
  
  
  
  
  
  
  cd ~/Downloads/
  
  load SC_AAL_90.mat
  
  areas = 90;
  C = SC/max(max(SC))*0.2;
  
  time_steps = 600;
  
  a = 0.2;
  we = 0.8;
  sig = 0.01;
  omega = 10;
  
  threshold = 0;
  polarity = 1;
  
  dt = 0.0025;
  
  xs = zeros(time_steps/2,areas);
  x = 0.1 * randn(areas,1);
  y = 0.1 * randn(areas,1);
  nn = 0;
  
  Isubdiag = find(tril(ones(areas),-1));
  
  
  for t=0:dt:500
    t
    sumax = C*x-diag(C*repmat(x',areas,1));
    sumay = C*y-diag(C*repmat(y',areas,1));
    
    x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
    y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
    
  end
  
  for t=0:dt:5000  %32000
    t
    % if t > 500
    %     a = 0.1;
    % end
    
    sumax = C*x-diag(C*repmat(x',areas,1));
    sumay = C*y-diag(C*repmat(y',areas,1));
    
    x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
    y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
    
    
    if mod(t,2)==0
      nn=nn+1;
      xs(nn,:)=x';
      ri=sqrt(x.*x+y.*y);
      rimax=max(ri);
      ri=ri/rimax;
      kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
    end
    
  end
  
  
  FC_simul = corrcoef(xs(1:nn,:));
  
  
  % RUN DFA
  f_win = [1 100];
  c_win = [0.5 150];
  
  siginfo = nbt_Info;
  siginfo.converted_sample_frequency =  5;
  
  t = nbt_doDFA(abs(hilbert(xs)), siginfo, [1 100],[0.5 150],0.5,0,0,[]);
  % cc = corrcoef(FC_emp(Isubdiag),FC_simul(Isubdiag));
  v = var(abs(hilbert(xs)),[],1);
  m = mean(abs(hilbert(xs)));
  
  
  %%
  % -0.2
  m = 0.0117
  v = 3.78e-05
  d = 0.56
  
  % 0
  m = 0.271
  v = 2.08e-04
  d = 0.92
  
  % 0.2
  m = 0.44;
  v = 8.72e-05;
  d = 0.56
  
  save('~/pmod/proc/pmod_hopf_par2.mat','t','v','m')
  
  
  Corr = cc(1,2);
  
  Corr_sim_mean = mean(mean(FC_simul));
  
  data = xs';
  [Peak,Amp] = Threshold_fMRI(data, areas, threshold, polarity);
  
  figure;
  subplot(2,1,1)
  pcolor(Peak)
  subplot(2,1,2)
  plot(1:time_steps/2+1, xs(:,1:66))
  xlim([0 time_steps/2+1])
  
