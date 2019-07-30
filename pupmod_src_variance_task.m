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
% addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

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

% t = license('test','signal_toolbox');
% if t
% %   continue
% else
%   error('Toolbox not available');
% end
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % %
    if ~exist(sprintf([outdir 'pupmod_src_variance_task_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_src_variance_task_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
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
      end
    end
    %
    save(sprintf([outdir 'pupmod_src_variance_task_s%d_m%d_v%d.mat'],isubj,m,v),'outp');
    
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
%% LOAD VARIANCE DATA

ord   = pconn_randomization;
for m = 1 : 3
  for isubj = SUBJLIST
    isubj
    for ifoi = 1:13
      for iblock = 1 : 2
        
        im = find(ord(isubj,:)==m);
        
        load(sprintf([outdir 'pupmod_src_variance_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v));
        
        var_all_rest(:,isubj,m,ifoi,iblock) = var; clear var
        
        load(sprintf([outdir 'pupmod_src_variance_task_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v));
        
        var_all_task(:,isubj,m,ifoi,iblock) = var; clear var
        
      end
    end
  end
end

var_all_rest = nanmean(var_all_rest(:,SUBJLIST,:,:,:),5);
var_all_task = nanmean(var_all_task(:,SUBJLIST,:,:,:),5);

%% EMPIRICAL

foi_range = [2 3 4 6 8 11 16 23 32 45 64 91 128];

for ifoi = 1 : 13
  
  [h,~,~,s]=ttest(var_all_rest(:,:,2,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_atx_pos_rest(ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_rest(ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_rest(:,:,3,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
  n_dpz_pos_rest(ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
  n_dpz_neg_rest(ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_task(:,:,2,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
  n_atx_pos_task(ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
  n_atx_neg_task(ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
  
  [h,~,~,s]=ttest(var_all_task(:,:,3,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
  n_dpz_pos_task(ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
  n_dpz_neg_task(ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
  
end

% PERMUTATION
nperm = 10000;
all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);

dat_cnt1 = var_all_task(:,:,[1 2],:);
dat_res1 = var_all_rest(:,:,[1 2],:);
dat_cnt2 = var_all_task(:,:,[1 3],:);
dat_res2 = var_all_rest(:,:,[1 3],:);

for iperm = 1 : nperm
  
  % within subjects permutation test
  fprintf('Perm #%d ...\n',iperm);
  
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    
    permdat_cnt1(:,i,1,:) = dat_cnt1(:,i,idx1(i),:);
    permdat_cnt1(:,i,2,:) = dat_cnt1(:,i,idx2(i),:);
    
    permdat_res1(:,i,1,:) = dat_res1(:,i,idx1(i),:);
    permdat_res1(:,i,2,:) = dat_res1(:,i,idx2(i),:);
    
  end
  
  for ifoi = 1 : 13
    
    h=ttest(permdat_res1(:,:,2,ifoi),permdat_res1(:,:,1,ifoi),'dim',2);
    n_atx_rest_perm(iperm,ifoi) = sum(h)./ length(h);
    
    h=ttest(permdat_cnt1(:,:,2,ifoi),permdat_cnt1(:,:,1,ifoi),'dim',2);
    n_atx_task_perm(iperm,ifoi) = sum(h)./ length(h);
    
  end
  
  for i = 1 : length(idx1)
    
    permdat_cnt1(:,i,1,:) = dat_cnt2(:,i,idx1(i),:);
    permdat_cnt1(:,i,2,:) = dat_cnt2(:,i,idx2(i),:);
    
    permdat_res1(:,i,1,:) = dat_res2(:,i,idx1(i),:);
    permdat_res1(:,i,2,:) = dat_res2(:,i,idx2(i),:);
    
  end
  
  for ifoi = 1 : 13
    
    [h,~,~,s]=ttest(var_all_rest(:,:,2,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
    perm.n_atx_pos_rest(iperm,ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_atx_neg_rest(iperm,ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(var_all_rest(:,:,3,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
    perm.n_dpz_pos_rest(iperm,ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_dpz_neg_rest(iperm,ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(var_all_task(:,:,2,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
    perm.n_atx_pos_task(iperm,ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_atx_neg_task(iperm,ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
    
    [h,~,~,s]=ttest(var_all_task(:,:,3,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
    perm.n_dpz_pos_task(iperm,ifoi) = sum((h>0)&(s.tstat>0))./ length(h);
    perm.n_dpz_neg_task(iperm,ifoi) = sum((h>0)&(s.tstat<0))./ length(h);
    
    
  end
end

for ifoi = 1 : 13
  
  p_atx_rest(ifoi) = 1-(sum(n_atx_rest(:,ifoi)>n_atx_rest_perm(:,ifoi))./nperm);
  p_atx_task(ifoi) = 1-(sum(n_atx_task(:,ifoi)>n_atx_task_perm(:,ifoi))./nperm);
  
  p_dpz_rest(ifoi) = 1-(sum(n_dpz_rest(:,ifoi)>n_dpz_rest_perm(:,ifoi))./nperm);
  p_dpz_task(ifoi) = 1-(sum(n_dpz_task(:,ifoi)>n_dpz_task_perm(:,ifoi))./nperm);
  
end

%%
figure; set(gcf,'color','w')

subplot(4,2,1); hold on
plot(n_atx_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_atx_neg_rest,'linewidth',2,'color',[0.2 0.5 1])
% plot(prctile(n_atx_rest_perm,95),'linewidth',1,'color',[1 0.5 0.2],'linestyle',':')
axis([0 14 -0.02 0.32]);
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
% xlabel('Frequency [Hz]');
ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,3); hold on
plot(n_atx_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_atx_neg_task,'linewidth',2,'color',[0.2 0.5 1])
% plot(prctile(n_atx_task_perm,95),'linewidth',1,'color',[1 0.1 0.1],'linestyle',':')
axis([0 14 -0.02 0.32]);
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
xlabel('Carrier frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

subplot(4,2,2); hold on
plot(n_dpz_pos_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_dpz_neg_rest,'linewidth',2,'color',[0.2 0.5 1])
% plot(prctile(n_dpz_rest_perm,95),'linewidth',1,'color',[0.2 0.5 1],'linestyle',':')

axis([0 14 -0.02 0.32]);
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
% xlabel('Frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Donepezil vs. placebo')
tp_editplots

subplot(4,2,4); hold on
plot(n_dpz_pos_task,'linewidth',2,'color',[1 0.5 0.2])
plot(n_dpz_neg_task,'linewidth',2,'color',[0.2 0.5 1])
% plot(prctile(n_dpz_task_perm,95),'linewidth',1,'color',[0.1 0.1 1],'linestyle',':')
axis([0 14 -0.02 0.32]);
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',num2cell([0 0.1 0.2 0.3]))
xlabel('Carrier frequency [Hz]'); %ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
% title('Atomoxetine vs. placebo')
tp_editplots

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_variance.eps'))








