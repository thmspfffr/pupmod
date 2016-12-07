function tc = pupmod_src_timecourse(dat,para)

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

% ------------------------------------------
% COMPUTE DFA IN SOURCE SPACE
% ------------------------------------------

if strcmp(para.filt,'lcmv')
  if strcmp(para.grid,'xcoarse')
    load('~/pconn/matlab/aalmask_grid_xcoarse.mat')
    load(para.sa);
    [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3));
  elseif strcmp(para.grid,'coarse')
    load('~/pconn/matlab/aalmask_grid_coarse.mat')
    load(para.sa);
    [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3));
  elseif strcmp(para.grid,'cortex')
    load(para.sa);
    [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3));
    A1 = getdipdir(nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3),A1);
  end
  
elseif strcmp(para.filt,'eloreta')
  if strcmp(para.grid,'xcoarse')
    load('~/pconn/matlab/aalmask_grid_xcoarse.mat')
    % v2
    load(para.sa);
    A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
    A1 = getdipdir(nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3),A1);
  elseif strcmp(para.grid,'coarse')
    load('~/pconn/matlab/aalmask_grid_coarse.mat')
    % v1
    load(para.sa);
    A1 = mkfilt_eloreta_v2(sa.L_coarse);
    A1 = getdipdir(nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3),A1);
  elseif strcmp(para.grid,'cortex')
    % v3
    load(para.sa);
    % watch out! L_coarse is correct even for cortex grid
    A1 = mkfilt_eloreta_v2(sa.L_coarse);
    A1 = getdipdir(nanmean(para.cs(:,:,para.foi(1,1):para.foi(1,2)),3),A1);
  end
end

siginfo = nbt_Info;
siginfo.converted_sample_frequency = para.fs;

% compute bp-filtered signal
ampenv = nbt_filter_fir(dat,para.foi(1,1),para.foi(1,2),siginfo.converted_sample_frequency,2/para.foi(1,1));

% project bp-filtered signal into source space
tc = ampenv*A1;

if para.aal
  
  for ireg = 1 : 91
    
    idx1           = aalgrid.mask==ireg;
    tc_aal(ireg,:) = nanmean(tc(:,idx1),2);
    
  end
  
  tc = tc_aal;
  
end





