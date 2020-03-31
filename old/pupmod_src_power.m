
%% pp_src_powcorr_test

% Call pupmod_all_powcorr_periphereal.m next, in order to clean
% estimated FCs from artifacts.

clear
% --------------------------------------------------------
% VERSION 23 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v               = 23;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = 2.^[1:.25:6];
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
fsample         = 400;
segleng         = 80;
segshift        = 40;
width           = 4;
weighting = 0;
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
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults()
outdir   = '/home/tpfeffer/pp/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

ord   = pconn_randomization;
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    
%     if ~exist(sprintf([outdir 'pupmod_src_power_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' outdir sprintf('pupmod_src_power_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end

    if ~weighting
      if v~=25
        powcorr = zeros(400,400,2,length(foi_range));
      else
        powcorr = zeros(90,90,2,length(foi_range));
      end
    else
      powcorr = zeros(90,90,2,length(foi_range));
    end
    for iblock = 1:2
      
      pars = [];
      pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      sa        = load(pars.sa);
      
      fprintf('Loading MEG data ...\n');
      
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
          if strcmp(allpara.grid,'cortex_lowres')
            powcorr(:,:,iblock,:) = nan(400,400,1,length(foi_range));
          else
            powcorr(:,:,iblock,:) = nan(90,90,1,length(foi_range));
          end
          continue
        else
          error('Data corrupt?')
        end
      end
      %       data.trial{1} = data.trial;
      %       dat = megdata2mydata(data);
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
        
        outp.pow(:,ifoi,iblock) = mean(abs(datasf1).^2,2);
        
        
%         fprintf('Computing pow corr  ...\n')
%         if ~weighting
%           powcorr(:,:,iblock,ifoi) = tp_data2orthopowcorr_wavelet(data,foi_range(ifoi),sa);
%         else
%           powcorr(:,:,iblock,ifoi) = tp_data2orthopowcorr_wavelet_weighted(data,foi_range(ifoi),sa);
%         end
%         
        
      end
           
    end
    save(sprintf([outdir 'pupmod_src_power_s%d_m%d_v%d.mat'],isubj,m,v),'outp');

  end
end
error('!')

%%
v=12;
clear allpowcorr
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      for ifoi = 1 :13
        try
          load(sprintf(['/home/tpfeffer/pp/proc/conn/' 'pp_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
          allpowcorr(:,:,isubj,m,iblock,ifoi,1) = single(squeeze(nanmean(powcorr,3))); clear powcorr
          %     load(sprintf([outdir 'pp_task_src_powcorr_test_s%d_m%d_v%d.mat'],isubj,im,v));
          %     allpowcorr(:,:,isubj,m,:,2) = single(squeeze(nanmean(powcorr,3))); clear powcorr
        catch me
          warning (sprintf('s%d m%d b%d',isubj,im,iblock))
          allpowcorr(:,:,isubj,m,iblock,ifoi,1) = nan(400,400); clear powcorr
          %   allpowcorr(:,:,isubj,m,:,2) = nan(400,400,3); clear powcorr
        end
      end
    end
  end
end

allpowcorr = squeeze(allpowcorr(:,:,SUBJLIST,:,:,:));

%% PLOT

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

close all

ifoi = 6:10; icond = 3; icont = 1;

% cmap =  hot
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);

par = nanmean(nanmean(nanmean(allpowcorr(:,:,:,icond,ifoi,icont),3),2),5)-nanmean(nanmean(nanmean(allpowcorr(:,:,:,1,ifoi,icont),3),2),5)%-nanmean(nanmean(p1(:,:,:,1,6),3),2);
[h] = ttest(squeeze(nanmean(nanmean(allpowcorr(:,:,:,icond,ifoi,icont),2),5)),squeeze(nanmean(nanmean(allpowcorr(:,:,:,1,ifoi,icont),2),5)),'dim',2,'alpha',0.025);%-nanmean(nanmean(p1(:,:,:,1,6),3),2)% par = par.*h;
% par(outp_atx.pval_p_atx(:,icond,ifoi)>=0.025) = 0;
par=  par.*h;
% cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-max([abs(min(par)) abs(max(par))]) max([abs(min(par)) abs(max(par))])];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)

%%
for ifoi = 6:10
  
  [h,~,~,s]=ttest(allpowcorr(:,:,:,3,ifoi),allpowcorr(:,:,:,1,ifoi),'dim',3)
  
  n_pos(ifoi) = 100*sum(sum((h>0)&(s.tstat>0)))/(400*400-400);
  n_neg(ifoi) = 100*sum(sum((h>0)&(s.tstat<0)))/(400*400-400);
  
end
%% PERMUTATION TEST


%   dat_res1 = single(squeeze(allpowcorr(:,:,:,[1 2],:)));
%   dat_res2 = single(squeeze(allpowcorr(:,:,:,[1 3],:)));
all_idx1 = randi(2,[size(SUBJLIST,2),1000]);

for iperm = 1 : 1000
  iperm
  % within subjects permutation test
  fprintf('Perm #%d\n',iperm);
  
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  %       for ifoi = 1 : 10
  %         ifoi
  for i = 1 : length(idx1)
    permdat_cnt1(:,:,i,1,:) = squeeze(allpowcorr(:,:,i,idx1(i),1:10));
    permdat_cnt1(:,:,i,2,:) = squeeze(allpowcorr(:,:,i,idx2(i),1:10));
  end
  
  [h,~,~,s]=ttest(permdat_cnt1(:,:,:,2,:),permdat_cnt1(:,:,:,1,:),'dim',3);
  
  n_pos_perm(iperm,:) = 100*sum(sum((h>0)&(s.tstat>0)))/(400*400-400);
  n_neg_perm(iperm,:) = 100*sum(sum((h>0)&(s.tstat<0)))/(400*400-400);
  
end
















