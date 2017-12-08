%% pupmod_src_amplitude_task

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
      if ~exist(sprintf([outdir 'pupmod_src_amplitude_task_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_amplitude_task_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
      for iblock = 1:1
        
        pars = [];
        pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        sa        = load(pars.sa);
        
        fprintf('Loading MEG data ...\n');

        try 
          load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          if ~exist(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
%             powcorr = nan(size(sa.sa.L_xcoarse,2),size(sa.sa.L_xcoarse,2));
            var = nan(90,90);
            save(sprintf([outdir 'pupmod_src_variance_task_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'var');
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

         [var] = tp_src_amp(dat,pars,sa);

        if size(var,1) < 100 && size(var,1) > 80
          pars = [];
          pars.grid = 'medium';
          var = tp_match_aal(pars,var);
        end
%    
       save(sprintf([outdir 'pupmod_src_amplitude_task_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'var');
        
      end
    end
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
%%
ord   = pconn_randomization;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        
        im = find(ord(isubj,:)==m); isubj

        load(sprintf([outdir 'pupmod_src_variance_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v));
        
        var_all_rest(:,isubj,m,ifoi,iblock) = diag(var); clear var
        
        load(sprintf([outdir 'pupmod_src_variance_task_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v));
        
        var_all_task(:,isubj,m,ifoi,iblock) = diag(var); clear var
        
      end
    end
  end
end

var_all_rest = nanmean(var_all_rest(:,SUBJLIST,:,:,:),5);
var_all_task = nanmean(var_all_task(:,SUBJLIST,:,:,:),5);

%%

% EMPIRICAL

foi_range = [2 3 4 6 8 11 16 23 32 45 64 91 128];

for ifoi = 1 : 13
  
   h=ttest(var_all_rest(:,:,2,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
   n_atx_rest(ifoi) = sum(h)./ length(h);
   
   h=ttest(var_all_rest(:,:,3,ifoi),var_all_rest(:,:,1,ifoi),'dim',2);
   n_dpz_rest(ifoi) = sum(h)./ length(h);
   
   h=ttest(var_all_task(:,:,2,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
   n_atx_task(ifoi) = sum(h)./ length(h);
   
   h=ttest(var_all_task(:,:,3,ifoi),var_all_task(:,:,1,ifoi),'dim',2);
   n_dpz_task(ifoi) = sum(h)./ length(h);
     
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
  
  disp(sprintf('Perm #%d',iperm));
  
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
    
    
    h=ttest(permdat_res1(:,:,2,ifoi),permdat_res1(:,:,1,ifoi),'dim',2);
    n_dpz_rest_perm(iperm,ifoi) = sum(h)./ length(h);
    
    h=ttest(permdat_cnt1(:,:,2,ifoi),permdat_cnt1(:,:,1,ifoi),'dim',2);
    n_dpz_task_perm(iperm,ifoi) = sum(h)./ length(h);
    
    
  end
  
  
  
end
  
for ifoi = 1 : 13
  
  p_atx_rest(ifoi) = 1-(sum(n_atx_rest(:,ifoi)>n_atx_rest_perm(:,ifoi))./nperm);
  p_atx_task(ifoi) = 1-(sum(n_atx_task(:,ifoi)>n_atx_task_perm(:,ifoi))./nperm);

  p_dpz_rest(ifoi) = 1-(sum(n_dpz_rest(:,ifoi)>n_dpz_rest_perm(:,ifoi))./nperm);
  p_dpz_task(ifoi) = 1-(sum(n_dpz_task(:,ifoi)>n_dpz_task_perm(:,ifoi))./nperm);


end

%%
figure; 

subplot(1,2,1); hold on
plot(n_atx_rest,'linewidth',2,'color',[1 0.5 0.2])
plot(n_atx_task,'linewidth',2,'color',[1 0.1 0.1])

plot(prctile(n_atx_rest_perm,95),'linewidth',1,'color',[1 0.5 0.2],'linestyle',':')
plot(prctile(n_atx_task_perm,95),'linewidth',1,'color',[1 0.1 0.1],'linestyle',':')


axis([0 14 -0.02 0.3]); axis square
set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
title('Atomoxetine vs. placebo')

subplot(1,2,2); hold on
plot(n_dpz_rest,'linewidth',2,'color',[0.2 0.5 1])
plot(n_dpz_task,'linewidth',2,'color',[0.1 0.1 1],'linestyle','-')

plot(prctile(n_dpz_rest_perm,95),'linewidth',1,'color',[0.2 0.5 1],'linestyle',':')
plot(prctile(n_dpz_task_perm,95),'linewidth',1,'color',[0.1 0.1 1],'linestyle',':')


axis([0 14 -0.02 0.3]); axis square
set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel(sprintf('Fraction of nodes \n with altered variance [%%]'))
title('Donepezil vs. placebo')

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_variance.pdf'))







  
