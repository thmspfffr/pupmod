%% pupmod_src_powcorr_permtest
% counts number of altered correlations and computes permututation test
% see hawellek et al. (2013) for details
% cleaned data is generated in pupmod_all_powcorr_periphereal

% data is collected in: pupmod_all_powcorr_plot_alteredcorr
% data is plotted in: pupmod_all_powcorr_plot_alteredcorr.m

% last update: 26-10-2018

clear
v = 1;
% alp: standard for all 0.05!!! except verison 13
nperm = 10000; alp = 0.05;
nfoi = 17;
par.subs = 250;
par.allperms = nperm/par.subs;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

cleandat = single(pupmod_loadpowcorr(v,SUBJLIST,1));

nvox = size(cleandat,1)*size(cleandat,1)-size(cleandat,1);

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
grid = sa.grid_cortex_lowres;

lab=tp_aal2yeo(grid);
%%

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')

if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_yeo_perms_subs%d_nperm%d_v%d.mat'],par(1).subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_yeo_subs%d_nperm%d_v%d.mat'],par(1).subs,nperm,v),'all_idx1');
else
  load(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_yeo_subs%d_nperm%d_v%d.mat'],par(1).subs,nperm,v));
end

for iperm = 1 : par(1).allperms
  
  if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_yeo_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' outdir sprintf('pupmod_src_powcorr_permtest_yeo_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par(1).subs
    
    iiperm = (iperm-1)*par(1).subs+kperm;
    
    % within subjects permutation test
    fprintf('Perm #%d\n',kperm);
    
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat_cnt1(:,:,i,1,:) = cleandat(:,:,i,idx1(i),2,:);
      permdat_cnt1(:,:,i,2,:) = cleandat(:,:,i,idx2(i),2,:);
      permdat_res1(:,:,i,1,:) = cleandat(:,:,i,idx1(i),1,:);
      permdat_res1(:,:,i,2,:) = cleandat(:,:,i,idx2(i),1,:);
    end
    
    % should this be the same as for ATX??
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    idx1(idx1==2)=3;
    idx2(idx2==2)=3;
    
    for i = 1 : length(idx1)
      permdat_cnt2(:,:,i,1,:) = cleandat(:,:,i,idx1(i),2,:);
      permdat_cnt2(:,:,i,2,:) = cleandat(:,:,i,idx2(i),2,:);
      permdat_res2(:,:,i,1,:) = cleandat(:,:,i,idx1(i),1,:);
      permdat_res2(:,:,i,2,:) = cleandat(:,:,i,idx2(i),1,:);
    end
    
    
    % -----------
    % compute ttest during task and atomoxetine
    [t_cnt1,~,~,s] = ttest(atanh(permdat_cnt1(:,:,:,2,:)),atanh(permdat_cnt1(:,:,:,1,:)),'dim',3,'alpha',alp);
    t_cnt1 = squeeze(t_cnt1.*sign(s.tstat)); clear s
    % compute ttest during rest and atomoxetine
    [t_res1,~,~,s] = ttest(atanh(permdat_res1(:,:,:,2,:)),atanh(permdat_res1(:,:,:,1,:)),'dim',3,'alpha',alp);
    t_res1 = squeeze(t_res1.*sign(s.tstat)); clear s
    % compute ttest during task and donepezil
    [t_cnt2,~,~,s] = ttest(atanh(permdat_cnt2(:,:,:,2,:)),atanh(permdat_cnt2(:,:,:,1,:)),'dim',3,'alpha',alp);
    t_cnt2 = squeeze(t_cnt2.*sign(s.tstat)); clear s
    % compute ttest during rest and donepezil
    [t_res2,~,~,s] = ttest(atanh(permdat_res2(:,:,:,2,:)),atanh(permdat_res2(:,:,:,1,:)),'dim',3,'alpha',alp);
    t_res2 = squeeze(t_res2.*sign(s.tstat)); clear s
    
    for ilab=1:7
      nvox = (sum(lab==ilab)*sum(lab==ilab))-sum(lab==ilab);
      % -----------------------
      % NUMBER OF ALTERED CORRELATONS - across space
      % -----------------------
      % count number of altered connections (atx, task)
      par(ilab).tperm_cnt1_n(kperm,:)=nansum(nansum(t_cnt1(lab==ilab,lab==ilab,:)<0))./nvox;
      par(ilab).tperm_cnt1_p(kperm,:)=nansum(nansum(t_cnt1(lab==ilab,lab==ilab,:)>0))./nvox;
      % count number of altered connections (atx, rest)
      par(ilab).tperm_res1_n(kperm,:)=nansum(nansum(t_res1(lab==ilab,lab==ilab,:)<0))./nvox;
      par(ilab).tperm_res1_p(kperm,:)=nansum(nansum(t_res1(lab==ilab,lab==ilab,:)>0))./nvox;
      % count number of altered connections (dpz, task)
      par(ilab).tperm_cnt2_n(kperm,:)=nansum(nansum(t_cnt2(lab==ilab,lab==ilab,:)<0))./nvox;
      par(ilab).tperm_cnt2_p(kperm,:)=nansum(nansum(t_cnt2(lab==ilab,lab==ilab,:)>0))./nvox;
      % count number of altered connections (dpz, rest)
      par(ilab).tperm_res2_n(kperm,:)=nansum(nansum(t_res2(lab==ilab,lab==ilab,:)<0))./nvox;
      par(ilab).tperm_res2_p(kperm,:)=nansum(nansum(t_res2(lab==ilab,lab==ilab,:)>0))./nvox;
      
    end
    
    
  end
  clear permdat_cnt1 permdat_cnt2 permdat_res1 permdat_res2
  
  
  
  save(sprintf([outdir 'pupmod_src_powcorr_permtest_yeo_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'par')
  
  try
    pause(randi(3))
    load(sprintf([outdir 'pupmod_src_powcorr_permtest_yeo_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v))
  catch me
    save(sprintf([outdir 'pupmod_src_powcorr_permtest_yeo_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'par')
  end
  
end
%%
