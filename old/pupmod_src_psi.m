%% pupmod_src_powcorr

clear 

% --------------------------------------------------------
% VERSION 1 - WEIGHTED AAL
% --------------------------------------------------------
v               = 1;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal_6mm';
foi_range       = 10;
para.segleng    = 1600;
para.bpfreq     = [8 12; 9 14];
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

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1 : length(foi_range)
% %         
      if ~exist(sprintf([outdir 'pupmod_src_psi_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_psi_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
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
        
        
        cs = data2cs_event(dat,1600,800,size(dat,1),49);
        f = 0:0.25:12; 
        cs = cs(:,:,find(f==8):find(f==12));
        
        psi = tp_psi(dat,pars,sa,cs);
       

        if size(psi,1) < 100 && size(psi,1) > 80
          pars = [];
          pars.grid = 'medium';
          psi = tp_match_aal(pars,psi);
        end
%    
       save(sprintf([outdir 'pupmod_src_psi_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'psi');
        
      end
    end
  end
end
  
error('!')

%%
ord   = pconn_randomization;

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1 : length(foi_range)
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
      
       load(sprintf([outdir 'pupmod_src_psi_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
       
       psi_all(:,:,isubj,m,ifoi,iblock) = psi;
       
      end
      
    end
  end
end


psi_all = squeeze(psi_all(:,:,SUBJLIST,:,:,:));

%% LOOK SPECIFICALLY AT THOSE CONNECTIONS THAT ARE DIFFERENT IN THE FIRST PLACE!

[h,~,~,s]=ttest(psi_all(:,:,:,2),psi_all(:,:,:,1),'dim',3);

n_p_atx = sum(sum(triu(h.*sign(s.tstat),1)>0))/((90*90-90)/2);
n_n_atx = sum(sum(triu(h.*sign(s.tstat),1)<0))/((90*90-90)/2);

[h,~,~,s]=ttest(psi_all(:,:,:,3),psi_all(:,:,:,1),'dim',3);

n_p_dpz = sum(sum(triu(h.*sign(s.tstat),1)>0))/((90*90-90)/2);
n_n_dpz = sum(sum(triu(h.*sign(s.tstat),1)<0))/((90*90-90)/2);

nperm = 2000;

for iperm = 1 : nperm
  
  iperm
  
  idx(:,1) = randi(2,[28 1]);
  idx(:,2) = 3-idx(:,1);
  
  for isubj = 1 : 28
    
    permdat(:,:,isubj,1) = psi_all(:,:,isubj,idx(isubj,1));
    permdat(:,:,isubj,2) = psi_all(:,:,isubj,idx(isubj,2));
    
  end
  
  [h,~,~,s]=ttest(permdat(:,:,:,2),permdat(:,:,:,1),'dim',3);

  
  n_p_atx_perm(iperm) = sum(sum(triu(h.*sign(s.tstat),1)>0))/((90*90-90)/2);
  n_n_atx_perm(iperm) = sum(sum(triu(h.*sign(s.tstat),1)<0))/((90*90-90)/2);

end

  
  


  











