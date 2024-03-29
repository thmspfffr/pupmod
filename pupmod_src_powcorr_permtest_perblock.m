%% pupmod_src_powcorr_permtest
% counts number of altered correlations and computes permututation test
% see hawellek et al. (2013) for details
% cleaned data is generated in pupmod_all_powcorr_periphereal

% data is collected in: pupmod_all_powcorr_plot_alteredcorr
% data is plotted in: pupmod_all_powcorr_plot_alteredcorr.m

% last update: 26-10-2018

clear
v = 12;
FOI = 6;
% alp: standard for all 0.05!!! except verison 13
nperm = 10000; alp = 0.05;

par.subs = 250;
par.allperms = nperm/par.subs;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
% addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

avg = 0;
cleandat = pupmod_loadpowcorr(v,avg);

nvox = size(cleandat,1)*size(cleandat,1)-size(cleandat,1);

%%
for iblock = 1 : 2
  
  tmp = clock;
  
  seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);
  
  rng(seed,'twister')
  
  if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_block%d_v%d.mat'],par.subs,nperm,iblock,v))
    all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
    save(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_block%d_v%d.mat'],par.subs,nperm,iblock,v),'all_idx1');
  else
    load(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_block%d_v%d.mat'],par.subs,nperm,iblock,v));
  end
  
  dat_res1 = single(squeeze(cleandat(:,:,:,[1 2],1,:,iblock)));
  dat_cnt1 = single(squeeze(cleandat(:,:,:,[1 2],2,:,iblock)));
  dat_cnt2 = single(squeeze(cleandat(:,:,:,[1 3],2,:,iblock)));
  dat_res2 = single(squeeze(cleandat(:,:,:,[1 3],1,:,iblock)));
  
  taskvsrest(:,:,:,1,:) = dat_res1(:,:,:,1,:);
  taskvsrest(:,:,:,2,:) = dat_cnt1(:,:,:,1,:);
  
  for iperm = 1 : par.allperms
    
    if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_block%d_v%d_processing.txt'],iperm,nperm,iblock,v))
      system(['touch ' outdir sprintf('pupmod_src_powcorr_permtest_iperm%d_block%d_nperm%d_v%d_processing.txt',iperm,nperm,iblock,v)]);
    else
      continue
    end
    
    for kperm = 1 : par.subs
      
      iiperm = (iperm-1)*par.subs+kperm;
      
      % within subjects permutation test
      fprintf('Perm #%d\n',kperm);
      
      idx1 = all_idx1(:,iiperm);
      idx2 = 3-idx1;
      
      for i = 1 : length(idx1)
        permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
        permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);
        permdat_res1(:,:,i,1,:) = dat_res1(:,:,i,idx1(i),:);
        permdat_res1(:,:,i,2,:) = dat_res1(:,:,i,idx2(i),:);
      end
      
      % should this be the same as for ATX??
      idx1 = all_idx1(:,iiperm);
      idx2 = 3-idx1;
      
      for i = 1 : length(idx1)
        permdat_cnt2(:,:,i,1,:) = dat_cnt2(:,:,i,idx1(i),:);
        permdat_cnt2(:,:,i,2,:) = dat_cnt2(:,:,i,idx2(i),:);
        permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
        permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:);
      end
      
      % TASK VS REST ------------------------------------------------------
      idx1 = all_idx1(:,iiperm);
      idx2 = 3-idx1;
      
      for i = 1 : length(idx1)
        taskvsrest_perm(:,:,i,1,:) = taskvsrest(:,:,i,idx1(i),:);
        taskvsrest_perm(:,:,i,2,:) = taskvsrest(:,:,i,idx2(i),:);
      end
      
      for ifoi = FOI
        
        % -----------
        % compute ttest during task and atomoxetine
        [t_cnt1,~,~,s] = ttest(atanh(permdat_cnt1(:,:,:,2,ifoi)),atanh(permdat_cnt1(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_cnt1 = t_cnt1.*sign(s.tstat); clear s
        % compute ttest during rest and atomoxetine
        [t_res1,~,~,s] = ttest(atanh(permdat_res1(:,:,:,2,ifoi)),atanh(permdat_res1(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_res1 = t_res1.*sign(s.tstat); clear s
        % compute ttest during task and donepezil
        [t_cnt2,~,~,s] = ttest(atanh(permdat_cnt2(:,:,:,2,ifoi)),atanh(permdat_cnt2(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_cnt2 = t_cnt2.*sign(s.tstat); clear s
        % compute ttest during rest and donepezil
        [t_res2,~,~,s] = ttest(atanh(permdat_res2(:,:,:,2,ifoi)),atanh(permdat_res2(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_res2 = t_res2.*sign(s.tstat); clear s
        
        % compute double contrast, test connections for context-dependence
        [t_all1,~,~,s] = ttest(atanh(permdat_res1(:,:,:,2,ifoi))-atanh(permdat_res1(:,:,:,1,ifoi)),atanh(permdat_cnt1(:,:,:,2,ifoi))-atanh(permdat_cnt1(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_all1 = t_all1.*sign(s.tstat); clear s
        % compute double contrast, test connections for context-dependence
        [t_all2,~,~,s] = ttest(atanh(permdat_res2(:,:,:,2,ifoi))-atanh(permdat_res2(:,:,:,1,ifoi)),atanh(permdat_cnt2(:,:,:,2,ifoi))-atanh(permdat_cnt2(:,:,:,1,ifoi)),'dim',3,'alpha',alp);
        t_all2 = t_all2.*sign(s.tstat); clear s
        % -----------------------
        % NUMBER OF ALTERED CORRELATONS - across space
        % -----------------------
        % count number of altered connections (atx, task)
        par.tperm_cnt1_n(kperm,ifoi)=nansum(nansum(t_cnt1<0))./nvox;
        par.tperm_cnt1_p(kperm,ifoi)=nansum(nansum(t_cnt1>0))./nvox;
        % count number of altered connections (atx, rest)
        par.tperm_res1_n(kperm,ifoi)=nansum(nansum(t_res1<0))./nvox;
        par.tperm_res1_p(kperm,ifoi)=nansum(nansum(t_res1>0))./nvox;
        % count number of altered connections (dpz, task)
        par.tperm_cnt2_n(kperm,ifoi)=nansum(nansum(t_cnt2<0))./nvox;
        par.tperm_cnt2_p(kperm,ifoi)=nansum(nansum(t_cnt2>0))./nvox;
        % count number of altered connections (dpz, rest)
        par.tperm_res2_n(kperm,ifoi)=nansum(nansum(t_res2<0))./nvox;
        par.tperm_res2_p(kperm,ifoi)=nansum(nansum(t_res2>0))./nvox;
        % number of altered connections, irrespective of direction (atx)
        par.tperm_atx_during_task(kperm,ifoi)=nansum(nansum(abs(t_cnt1)))./nvox;
        par.tperm_atx_during_rest(kperm,ifoi)=nansum(nansum(abs(t_res1)))./nvox;
        % number of altered connections, irrespective of direction (dpz)
        par.tperm_dpz_during_task(kperm,ifoi)=nansum(nansum(abs(t_cnt2)))./nvox;
        par.tperm_dpz_during_rest(kperm,ifoi)=nansum(nansum(abs(t_res2)))./nvox;
        
        % -----------------------
        % NUMBER OF ALTERED CORRELATONS - per voxel
        % -----------------------
        % count number of altered connections (atx, task)
        par.tperm_cnt1_pervoxel_n(:,kperm,ifoi)=nansum(t_cnt1<0)./size(cleandat,1);
        par.tperm_cnt1_pervoxel_p(:,kperm,ifoi)=nansum(t_cnt1>0)./size(cleandat,1);
        % count number of altered connections (atx, rest)
        par.tperm_res1_pervoxel_n(:,kperm,ifoi)=nansum(t_res1<0)./size(cleandat,1);
        par.tperm_res1_pervoxel_p(:,kperm,ifoi)=nansum(t_res1>0)./size(cleandat,1);
        % count number of altered connections (dpz, task)
        par.tperm_cnt2_pervoxel_n(:,kperm,ifoi)=nansum(t_cnt2<0)./size(cleandat,1);
        par.tperm_cnt2_pervoxel_p(:,kperm,ifoi)=nansum(t_cnt2>0)./size(cleandat,1);
        % count number of altered connections (dpz, rest)
        par.tperm_res2_pervoxel_n(:,kperm,ifoi)=nansum(t_res2<0)./size(cleandat,1);
        par.tperm_res2_pervoxel_p(:,kperm,ifoi)=nansum(t_res2>0)./size(cleandat,1);
        % number of altered connections, irrespective of direction (atx)
        par.tperm_atx_during_task_pervoxel(:,kperm,ifoi)=nansum(abs(t_cnt1))./size(cleandat,1);
        par.tperm_atx_during_rest_pervoxel(:,kperm,ifoi)=nansum(abs(t_res1))./size(cleandat,1);
        % number of altered connections, irrespective of direction (dpz)
        par.tperm_dpz_during_task_pervoxel(:,kperm,ifoi)=nansum(abs(t_cnt2))./size(cleandat,1);
        par.tperm_dpz_during_rest_pervoxel(:,kperm,ifoi)=nansum(abs(t_res2))./size(cleandat,1);
        
        % -----------------------
        % CONTEXT DEPENDENCE - across space
        % -----------------------
        % context dependence: diff in counts (atx) between task and rest
        par.tperm_context_diff_atx_n(kperm,ifoi)= par.tperm_cnt1_n(kperm,ifoi) - par.tperm_res1_n(kperm,ifoi);
        par.tperm_context_diff_atx_p(kperm,ifoi)= par.tperm_cnt1_p(kperm,ifoi) - par.tperm_res1_p(kperm,ifoi);
        % context dependence: diff in counts (dpz) between task and rest
        par.tperm_context_diff_dpz_n(kperm,ifoi)= par.tperm_cnt2_n(kperm,ifoi) - par.tperm_res2_n(kperm,ifoi);
        par.tperm_context_diff_dpz_p(kperm,ifoi)= par.tperm_cnt2_p(kperm,ifoi) - par.tperm_res2_p(kperm,ifoi);
        % context-dependence: diff in counts, irrespective of direction
        par.tperm_atx_context(kperm,ifoi)=par.tperm_atx_during_task(kperm,ifoi)-par.tperm_atx_during_rest(kperm,ifoi);
        par.tperm_dpz_context(kperm,ifoi)=par.tperm_dpz_during_task(kperm,ifoi)-par.tperm_dpz_during_rest(kperm,ifoi);
        
        % -----------------------
        % CONTEXT DEPENDENCE - per voxel
        % -----------------------
        % context dependence: diff in counts (atx) between task and rest
        par.tperm_context_diff_atx_n_pervoxel(:,kperm,ifoi) = par.tperm_cnt1_pervoxel_n(:,kperm,ifoi) - par.tperm_res1_pervoxel_n(:,kperm,ifoi);
        par.tperm_context_diff_atx_p_pervoxel(:,kperm,ifoi) = par.tperm_cnt1_pervoxel_p(:,kperm,ifoi) - par.tperm_res1_pervoxel_p(:,kperm,ifoi);
        % context dependence: diff in counts (dpz) between task and rest
        par.tperm_context_diff_dpz_n_pervoxel(:,kperm,ifoi) = par.tperm_cnt2_pervoxel_n(:,kperm,ifoi) - par.tperm_res2_pervoxel_n(:,kperm,ifoi);
        par.tperm_context_diff_dpz_p_pervoxel(:,kperm,ifoi) = par.tperm_cnt2_pervoxel_p(:,kperm,ifoi) - par.tperm_res2_pervoxel_p(:,kperm,ifoi);
        
        % test double dissociation (atx vs. dpz, rest vs. task)
        a = par.tperm_context_diff_atx_n(kperm,ifoi)+par.tperm_context_diff_atx_p(kperm,ifoi);
        b = par.tperm_context_diff_dpz_n(kperm,ifoi)+par.tperm_context_diff_dpz_p(kperm,ifoi);
        par.tperm_doubledissociation(kperm,ifoi) = a-b;
        
        % context-dependence: number of context-dependent connections
        par.tperm_atx_context_test(kperm,ifoi) = nansum(nansum(t_all1))./nvox;
        par.tperm_dpz_context_test(kperm,ifoi) = nansum(nansum(t_all2))./nvox;
        
        % TASK VS REST ------------------------------------------------------
        [t_tvsr,~,~,s] = ttest(taskvsrest_perm(:,:,:,2,ifoi),taskvsrest_perm(:,:,:,1,ifoi),'dim',3,'alpha',alp);
        t_tvsr = t_tvsr.*sign(s.tstat); clear s
        
        par.tperm_taskvsrest_p(kperm,ifoi) = nansum(nansum(t_tvsr>0))./nvox;
        par.tperm_taskvsrest_n(kperm,ifoi) = nansum(nansum(t_tvsr<0))./nvox;
        
        par.tperm_taskvsrest_pervox_p(:,kperm,ifoi) = nansum(t_tvsr>0)./size(cleandat,1);
        par.tperm_taskvsrest_pervox_n(:,kperm,ifoi) = nansum(t_tvsr<0)./size(cleandat,1);
        
      end
    end
    
    save(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_block%d_v%d.mat'],iperm,nperm,iblock,v),'par')
    
    try
      pause(randi(3))
      load(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_block%d_v%d.mat'],iperm,nperm,iblock,v))
    catch me
      save(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_block%d_v%d.mat'],iperm,nperm,iblock,v),'par')
    end
    
  end
end
%%
