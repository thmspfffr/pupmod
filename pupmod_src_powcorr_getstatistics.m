function [outp] = pupmod_src_powcorr_getstatistics(para)
%% OBTAIN CORRECTED STATISTICAL THRESHOLDS FROM PERMUTATION DISTRIBUTIONS
% This function obtains corrected statstical thresholds based on Hawellek
% et al. (2013) Journal of Neuroscience.
% After this function is executed, the permutation test can be run a second
% time in order to obtain permutation results of pooled, significant
% frequencies on single voxel level (when para.type == 'local').
% -----------------
% Takes the following inputs:
% -----------------
% para.type = 'global' / 'local'
% para.cond = 'atx' / 'dpz' / 'taskvsrest'
% para.ver = version, see pupmod_src_powcorr.m / pupmod_src_powcorr_permtest.m
% para.nperm = number of permutations (see pupmod_src_powcorr_permtest.m)
% para.nsubs = number of permutations (see pupmod_src_powcorr_permtest.m)
% para.alpha = global alpha level (default = 0.05)
% para.alp   = alpha level to get number of altered corr. (default =  0.05)
% para.correction_method = 'ranks' or 'single_threshold'
% para.empirical = empirical data (altered corr)
% ------------------
% GLOBAL STATISTICS - across all voxels (para.type = 'global')
% ------------------
%   (1a) Group labels are randomly permuted *nperm* times
%   (2a) Number of altered correlations is computed for each frequency and
%        effect direction (increases/decreases), resulting in a permutation
%        distribution of altered correlations D1.
%   (3a) For each frequency, this distribution is turned into ranks, resulting
%        in a permutation rankdistribution R1. From this, for each resample,
%        the maximum rank is determined *across* frequencies, resulting in a
%        maximum rank distribution R1_max.  T
%   (4a) This distribution is then used in order to derive corrected p-values,
%        for each frequency, from D1.
% ------------------
% LOCAL STATISTICS - per voxel (para.type = 'local')
% ------------------
%   (1b) Same as 1a (in fact, based on the same permutations)
%   (2b) Number of altered correlations is computed *per voxel* for each
%        frequency and effect direction, resulting in D2
%   (2c)
%
% ------------------
% Last edited 27-07/2018
% ------------------
% Things to add:
% - Double dissociation per voxel
%
% tpfeffer | thmspfffr@gmail.com | 2018

%% FIRST LOAD EMPIRICAL DATA AND COMPUTE EMPIRICAL EFFECTS
% Number of altered correlations across all voxels

allperms = para.nperm/para.nsubs;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

if ~isfield(para,'emp')
  fprintf('Loading data...\n')
  cleandat = pupmod_loadpowcorr(para.ver,SUBJLIST,1);
  emp = pupmod_compute_altered_correlations(cleandat,para);
  para.fcsize = size(cleandat,1);
else
  emp = para.emp;
  para.fcsize = size(emp.n_p_atx_pervoxel,1);
end

clear cleandat s_fc

% --------------------------
% LOAD PERMUTATION DISTRIBUTION
% --------------------------

if strcmp(para.type, 'global')
  for iperm = 1 : allperms
    
    fprintf('Load permutation distributions: %d / %d ...\n',iperm,allperms)
    
    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,para.nperm,para.ver),'par')
    
    % --------------
    % ATOMOXETINE: global effects
    % --------------
    % rest
    perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)   = par.tperm_res1_p;
    perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)   = par.tperm_res1_n;
    % task
    perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)   = par.tperm_cnt1_p;
    perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)   = par.tperm_cnt1_n;
    % irrespective of context
    perm_n_all_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1) = par.tperm_atx_during_rest;
    perm_n_all_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2) = par.tperm_atx_during_task;
    % context-dependence
    context_allconn_atx_p((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_p-par.tperm_cnt1_p;
    context_allconn_atx_n((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_n-par.tperm_cnt1_n;
    
    % --------------
    % DONEPEZIL
    % --------------
    % rest
    perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)   = par.tperm_res2_p;
    perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)   = par.tperm_cnt2_p;
    % task
    perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)   = par.tperm_res2_n;
    perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)   = par.tperm_cnt2_n;
    % irrespective of context
    perm_n_all_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1) = par.tperm_dpz_during_rest;
    perm_n_all_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2) = par.tperm_dpz_during_task;
    % context-dependence
    context_allconn_dpz_p((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_p-par.tperm_cnt2_p;
    context_allconn_dpz_n((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_n-par.tperm_cnt2_n;
    
    % TASK VS REST
    % --------------
    % global effects
    perm_taskvsrest_n_p((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_p;
    perm_taskvsrest_n_n((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_n;
    
    % double dissociation
    % --------------
    perm_doubledissociation((iperm-1)*par.subs+1:(iperm)*par.subs,:) = context_allconn_atx_p((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)-context_allconn_dpz_n((iperm-1)*par.subs+1:(iperm)*par.subs,:,1);
  end
  
  emp.doubledissociation_emp = emp.n_p_context_atx-emp.n_n_context_dpz;
  
  % ------------------
  % GET PERCENTILES FROM PERMUTATION DISTRIBUTION (use as shaded error bar)
  % ------------------
  % Atomoxetine
  outp.percentiles_atx_pos(:,:,1)  = [prctile(perm_n_p_atx(:,:,1),0.5); prctile(perm_n_p_atx(:,:,1),99.5)];
  outp.percentiles_atx_pos(:,:,2)  = [prctile(perm_n_p_atx(:,:,2),0.5); prctile(perm_n_p_atx(:,:,2),99.5)];
  outp.percentiles_atx_neg(:,:,1)  = [prctile(perm_n_n_atx(:,:,1),0.5); prctile(perm_n_n_atx(:,:,1),99.5)];
  outp.percentiles_atx_neg(:,:,2)  = [prctile(perm_n_n_atx(:,:,2),0.5); prctile(perm_n_n_atx(:,:,2),99.5)];
  outp.percentiles_atx_context_pos = [prctile(context_allconn_atx_p,0.5); prctile(context_allconn_atx_p,99.5)];
  outp.percentiles_atx_context_neg = [prctile(context_allconn_atx_n,0.5); prctile(context_allconn_atx_n,99.5)];
  % Donepezil
  outp.percentiles_dpz_pos(:,:,1)  = [prctile(perm_n_p_dpz(:,:,1),0.5); prctile(perm_n_p_dpz(:,:,1),99.5)];
  outp.percentiles_dpz_pos(:,:,2)  = [prctile(perm_n_p_dpz(:,:,2),0.5); prctile(perm_n_p_dpz(:,:,2),99.5)];
  outp.percentiles_dpz_neg(:,:,1)  = [prctile(perm_n_n_dpz(:,:,1),0.5); prctile(perm_n_n_dpz(:,:,1),99.5)];
  outp.percentiles_dpz_neg(:,:,2)  = [prctile(perm_n_n_dpz(:,:,2),0.5); prctile(perm_n_n_dpz(:,:,2),99.5)];
  outp.percentiles_atx_context_pos = [prctile(context_allconn_dpz_p,0.5); prctile(context_allconn_dpz_p,99.5)];
  outp.percentiles_atx_context_neg = [prctile(context_allconn_dpz_n,0.5); prctile(context_allconn_dpz_n,99.5)];
  % Double dissociation
  outp.percentiles_doubledissociaton = [prctile(perm_doubledissociation,0.5); prctile(perm_doubledissociation,99.5)];
  % ------------------
  
  %%
  fprintf('Obtaining corrected p-values ...\n')
      
  % -----------------------------------------------------------------------
  % SINGLE THRESHOLD PERMUTATION TEST
  % -----------------------------------------------------------------------
  if strcmp(para.correction_method,'single_threshold')
    
    % based on maximum statistic, i.e. max. fraction of altered correlations
    % across all frequencies. note that this could be biased against higher
    % frequencies which tend to be more local in nature.
    % see nichols & holmes (2001) hum brain mapp
    
    % ------------------
    % ATOMOXETINE
    % ------------------
    idx_R_res1_p      = max(abs(perm_n_p_atx(:,:,1)),[],2);
    idx_R_res1_n      = max(abs(perm_n_n_atx(:,:,1)),[],2);
    idx_R_cnt1_p      = max(abs(perm_n_p_atx(:,:,2)),[],2);
    idx_R_cnt1_n      = max(abs(perm_n_n_atx(:,:,2)),[],2);
    idx_R_res1_all    = max(abs(perm_n_all_atx(:,:,1)),[],2);
    idx_R_cnt1_all    = max(abs(perm_n_all_atx(:,:,2)),[],2);
    idx_R_context1_p  = max(abs(context_allconn_atx_p),[],2);
    idx_R_context1_n  = max(abs(context_allconn_atx_n),[],2);
    % ------------------
    % DONEPEZIL
    % ------------------
    idx_R_res2_p      = max(abs(perm_n_p_dpz(:,:,1)),[],2);
    idx_R_res2_n      = max(abs(perm_n_n_dpz(:,:,1)),[],2);
    idx_R_cnt2_p      = max(abs(perm_n_p_dpz(:,:,2)),[],2);
    idx_R_cnt2_n      = max(abs(perm_n_n_dpz(:,:,2)),[],2);
    idx_R_res2_all    = max(abs(perm_n_all_dpz(:,:,1)),[],2);
    idx_R_cnt2_all    = max(abs(perm_n_all_dpz(:,:,2)),[],2);
    idx_R_context2_p  = max(abs(context_allconn_dpz_p),[],2);
    idx_R_context2_n  = max(abs(context_allconn_dpz_n),[],2);
    % ------------------
    % TASK VS REST
    % ------------------
    idx_R_taskvsrest_p = max(abs(perm_taskvsrest_n_p),[],2);
    idx_R_taskvsrest_n = max(abs(perm_taskvsrest_n_n),[],2);
    % ------------------
    % DOUBLE DISSOCIATION
    % ------------------ 
    idx_R_doublediss = max(abs(perm_doubledissociation),[],2);
    % ------------------
    for ifreq = para.nfreq
      % ------------------
      % ATOMOXETINE
      % ------------------
      outp.p_res1_p(ifreq)        = 1-sum(abs(emp.n_p_atx(ifreq,1))>abs(idx_R_res1_p))/para.nperm;
      outp.p_res1_n(ifreq)        = 1-sum(abs(emp.n_n_atx(ifreq,1))>abs(idx_R_res1_n))/para.nperm;
      outp.p_cnt1_p(ifreq)        = 1-sum(abs(emp.n_p_atx(ifreq,2))>abs(idx_R_cnt1_p))/para.nperm;
      outp.p_cnt1_n(ifreq)        = 1-sum(abs(emp.n_n_atx(ifreq,2))>abs(idx_R_cnt1_n))/para.nperm;
      outp.p_cnt_atx_all(ifreq)   = 1-sum(abs(emp.atx(ifreq,2))>abs(idx_R_cnt1_all))/para.nperm;
      outp.p_res_atx_all(ifreq)   = 1-sum(abs(emp.atx(ifreq,1))>abs(idx_R_res1_all))/para.nperm;
      outp.p_context_atx_p(ifreq) = 1-sum(abs(emp.n_p_context_atx(ifreq))>abs(idx_R_context1_p))/para.nperm;
      outp.p_context_atx_n(ifreq) = 1-sum(abs(emp.n_n_context_atx(ifreq))>abs(idx_R_context1_n))/para.nperm;
      % ------------------
      % DONEPEZIL
      % ------------------
      outp.p_res2_p(ifreq)        = 1-sum(abs(emp.n_p_dpz(ifreq,1))>abs(idx_R_res2_p))/para.nperm;
      outp.p_res2_n(ifreq)        = 1-sum(abs(emp.n_n_dpz(ifreq,1))>abs(idx_R_res2_n))/para.nperm;
      outp.p_cnt2_p(ifreq)        = 1-sum(abs(emp.n_p_dpz(ifreq,2))>abs(idx_R_cnt2_p))/para.nperm;
      outp.p_cnt2_n(ifreq)        = 1-sum(abs(emp.n_n_dpz(ifreq,2))>abs(idx_R_cnt2_n))/para.nperm;
      outp.p_cnt_dpz_all(ifreq)   = 1-sum(abs(emp.dpz(ifreq,2))>idx_R_cnt2_all)/para.nperm;
      outp.p_res_dpz_all(ifreq)   = 1-sum(abs(emp.dpz(ifreq,1))>idx_R_res2_all)/para.nperm;
      outp.p_context_dpz_p(ifreq) = 1-sum(abs(emp.n_p_context_dpz(ifreq))>abs(idx_R_context2_p))/para.nperm;
      outp.p_context_dpz_n(ifreq) = 1-sum(abs(emp.n_n_context_dpz(ifreq))>abs(idx_R_context2_n))/para.nperm;
      % ------------------
      % TASK VS REST
      % ------------------
      outp.pval_taskvsrest_p_corr(ifreq) = 1-sum(abs(emp.taskvsrest_p(ifreq))>idx_R_taskvsrest_p)/para.nperm;
      outp.pval_taskvsrest_n_corr(ifreq) = 1-sum(abs(emp.taskvsrest_n(ifreq))>idx_R_taskvsrest_n)/para.nperm;
      % ------------------
      % DOUBLE DISSOCIATION
      % ------------------ 
      outp.pval_doublediss(ifreq) = 1-sum(abs(emp.doubledissociation_emp(ifreq))>idx_R_doublediss)/para.nperm;
         
    end
  end
  
  %% --------------------------
  % LOCAL ANALYSIS (VOXEL LEVEL)
  % --------------------------
  % Correction is done using a single threshold permutation procedure
  
elseif strcmp(para.type,'local')
  
  perm_n_p_atx_pervoxel           = zeros(para.nperm,para.fcsize,length(para.nfreq),2,'single');
  perm_n_n_atx_pervoxel           = zeros(para.nperm,para.fcsize,length(para.nfreq),2,'single');
  perm_n_p_context_atx_pervoxel   = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  perm_n_n_context_atx_pervoxel   = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  perm_n_p_dpz_pervoxel           = zeros(para.nperm,para.fcsize,length(para.nfreq),2,'single');
  perm_n_n_dpz_pervoxel           = zeros(para.nperm,para.fcsize,length(para.nfreq),2,'single');
  perm_n_p_context_dpz_pervoxel   = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  perm_n_n_context_dpz_pervoxel   = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  perm_taskvsrest_n_p_pervoxel    = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  perm_taskvsrest_n_n_pervoxel    = zeros(para.nperm,para.fcsize,length(para.nfreq),1,'single');
  
  for iperm = 1 : allperms
    
    fprintf('Concatenating permutation distribution: perm%d ...\n',iperm)
    
    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,para.nperm,para.ver),'par')
    
    % --------------
    % ATOMOXETINE: local effects
    % --------------
    % rest (indicated by *res1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_n,[2 1 3]);
    % task (indicated by *cnt1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_n,[2 1 3]);
    perm_n_p_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_p_pervoxel,[2 1 3]);
    perm_n_n_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_n_pervoxel,[2 1 3]);
    % --------------
    % DONEPEZIL: local effects
    % --------------
    % rest (indicated by *res2* ending)
    perm_n_p_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1)   = permute(par.tperm_res2_pervoxel_p,[2 1 3]);
    perm_n_n_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1)   = permute(par.tperm_res2_pervoxel_n,[2 1 3]);
    % task (indicated by *cnt2* ending)
    perm_n_p_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2)   = permute(par.tperm_cnt2_pervoxel_p,[2 1 3]);
    perm_n_n_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2)   = permute(par.tperm_cnt2_pervoxel_n,[2 1 3]);
    % context-dependence
    perm_n_p_context_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_dpz_p_pervoxel,[2 1 3]);
    perm_n_n_context_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_dpz_n_pervoxel,[2 1 3]);
    % --------------
    % Task vs Rest: local effects
    % --------------
    perm_taskvsrest_n_p_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_p,[2 1 3]);
    perm_taskvsrest_n_n_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_n,[2 1 3]);
    
  end

  
  if strcmp(para.correction_method, 'ranks')
    error('not implemented')
  
  elseif strcmp(para.correction_method, 'single_threshold')
    
%     THIS IS IN ORDER TO CORRECT ACROSS SPACE *AND* FREQUENCIES
%     --------------------------------------
%     in order to work, same block within loop has to be commented out
%     -------

    d_max_p_atx     = squeeze(max(reshape(perm_n_p_atx_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 2]),[],2));
    d_max_n_atx     = squeeze(max(reshape(perm_n_n_atx_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 2]),[],2));
    d_max_ctx_p_atx = squeeze(max(reshape(perm_n_p_context_atx_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));
    d_max_ctx_n_atx = squeeze(max(reshape(perm_n_n_context_atx_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));
          
    d_max_p_dpz     = squeeze(max(reshape(perm_n_p_dpz_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 2]),[],2));
    d_max_n_dpz     = squeeze(max(reshape(perm_n_n_dpz_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 2]),[],2));
    d_max_ctx_p_dpz = squeeze(max(reshape(perm_n_p_context_dpz_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));
    d_max_ctx_n_dpz = squeeze(max(reshape(perm_n_n_context_dpz_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));
          
    d_max_p_tvr  =squeeze(max(reshape(perm_taskvsrest_n_p_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));
    d_max_n_tvr  = squeeze(max(reshape(perm_taskvsrest_n_n_pervoxel,[para.nperm para.fcsize*length(para.nfreq) 1]),[],2));

%     --------------------------------------

    for ifreq = para.nfreq
%      
      for ivox = 1 : para.fcsize
        outp.pval_vox_p_atx(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_p_atx_pervoxel(ivox,ifreq,:))' > d_max_p_atx)./para.nperm;
        outp.pval_vox_n_atx(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_n_atx_pervoxel(ivox,ifreq,:))' > d_max_n_atx)./para.nperm;
        outp.pval_vox_p_atx_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_p_context_atx_pervoxel(ivox,ifreq))' > d_max_ctx_p_atx)./para.nperm;
        outp.pval_vox_n_atx_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_n_context_atx_pervoxel(ivox,ifreq))' > d_max_ctx_n_atx)./para.nperm;
      
        outp.pval_vox_p_dpz(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_p_dpz_pervoxel(ivox,ifreq,:))' > d_max_p_dpz)./para.nperm;
        outp.pval_vox_n_dpz(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_n_dpz_pervoxel(ivox,ifreq,:))' > d_max_n_dpz)./para.nperm;
        outp.pval_vox_p_dpz_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_p_context_dpz_pervoxel(ivox,ifreq))' > d_max_ctx_p_dpz)./para.nperm;
        outp.pval_vox_n_dpz_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_n_context_dpz_pervoxel(ivox,ifreq))' > d_max_ctx_n_dpz)./para.nperm;
      
        outp.pval_vox_p_tvr(ivox,:,ifreq)   = 1-sum(squeeze(emp.taskvsrest_p_pervoxel(ivox,ifreq,:))' > d_max_p_tvr)./para.nperm;
        outp.pval_vox_n_tvr(ivox,:,ifreq)   = 1-sum(squeeze(emp.taskvsrest_n_pervoxel(ivox,ifreq,:))' > d_max_n_tvr)./para.nperm;
      end  
    end  
  end
  
elseif strcmp(para.type,'second_level_local')
  
  perm_n_p_atx_pervoxel           = zeros(para.nperm,para.fcsize,17,2,'single');
  perm_n_n_atx_pervoxel           = zeros(para.nperm,para.fcsize,17,2,'single');
  perm_n_p_context_atx_pervoxel   = zeros(para.nperm,para.fcsize,17,1,'single');
  perm_n_n_context_atx_pervoxel   = zeros(para.nperm,para.fcsize,17,1,'single');
  perm_n_p_dpz_pervoxel           = zeros(para.nperm,para.fcsize,17,2,'single');
  perm_n_n_dpz_pervoxel           = zeros(para.nperm,para.fcsize,17,2,'single');
  perm_n_p_context_dpz_pervoxel   = zeros(para.nperm,para.fcsize,17,1,'single');
  perm_n_n_context_dpz_pervoxel   = zeros(para.nperm,para.fcsize,17,1,'single');
  perm_taskvsrest_n_p_pervoxel    = zeros(para.nperm,para.fcsize,17,1,'single');
  perm_taskvsrest_n_n_pervoxel    = zeros(para.nperm,para.fcsize,17,1,'single');
  
  for iperm = 1 : allperms
  
  fprintf('Concatenating permutation distribution: perm%d ...\n',iperm)
    
    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,para.nperm,para.ver),'par')
    
    % --------------
    % ATOMOXETINE: local effects
    % --------------
    % rest (indicated by *res1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_n,[2 1 3]);
    % task (indicated by *cnt1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_n,[2 1 3]);
    perm_n_p_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_p_pervoxel,[2 1 3]);
    perm_n_n_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_n_pervoxel,[2 1 3]);
    % --------------
    % DONEPEZIL: local effects
    % --------------
    % rest (indicated by *res2* ending)
    perm_n_p_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1)   = permute(par.tperm_res2_pervoxel_p,[2 1 3]);
    perm_n_n_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1)   = permute(par.tperm_res2_pervoxel_n,[2 1 3]);
    % task (indicated by *cnt2* ending)
    perm_n_p_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2)   = permute(par.tperm_cnt2_pervoxel_p,[2 1 3]);
    perm_n_n_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2)   = permute(par.tperm_cnt2_pervoxel_n,[2 1 3]);
    % context-dependence
    perm_n_p_context_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_dpz_p_pervoxel,[2 1 3]);
    perm_n_n_context_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_dpz_n_pervoxel,[2 1 3]);
    % --------------
    % Task vs Rest: local effects
    % --------------
    perm_taskvsrest_n_p_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_p,[2 1 3]);
    perm_taskvsrest_n_n_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_n,[2 1 3]);
    
  end
    
  if strcmp(para.correction_method, 'single_threshold')
    
    nfreq = 1:17;
    d_max_p_atx     = squeeze(max(reshape(perm_n_p_atx_pervoxel,[para.nperm para.fcsize*length(nfreq) 2]),[],2));
    d_max_n_atx     = squeeze(max(reshape(perm_n_n_atx_pervoxel,[para.nperm para.fcsize*length(nfreq) 2]),[],2));
    d_max_ctx_p_atx = squeeze(max(reshape(perm_n_p_context_atx_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));
    d_max_ctx_n_atx = squeeze(max(reshape(perm_n_n_context_atx_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));
    
    d_max_p_dpz     = squeeze(max(reshape(perm_n_p_dpz_pervoxel,[para.nperm para.fcsize*length(nfreq) 2]),[],2));
    d_max_n_dpz     = squeeze(max(reshape(perm_n_n_dpz_pervoxel,[para.nperm para.fcsize*length(nfreq) 2]),[],2));
    d_max_ctx_p_dpz = squeeze(max(reshape(perm_n_p_context_dpz_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));
    d_max_ctx_n_dpz = squeeze(max(reshape(perm_n_n_context_dpz_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));
    
    d_max_p_tvr  =squeeze(max(reshape(perm_taskvsrest_n_p_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));
    d_max_n_tvr  = squeeze(max(reshape(perm_taskvsrest_n_n_pervoxel,[para.nperm para.fcsize*length(nfreq) 1]),[],2));

    
    for ifreq = para.nfreq
   
      for ivox = 1 : para.fcsize
        outp.pval_vox_p_atx(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_p_atx_pervoxel(ivox,ifreq,:))' > d_max_p_atx)./para.nperm;
        outp.pval_vox_n_atx(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_n_atx_pervoxel(ivox,ifreq,:))' > d_max_n_atx)./para.nperm;
        outp.pval_vox_p_atx_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_p_context_atx_pervoxel(ivox,ifreq))' > d_max_ctx_p_atx)./para.nperm;
        outp.pval_vox_n_atx_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_n_context_atx_pervoxel(ivox,ifreq))' > d_max_ctx_n_atx)./para.nperm;
      
        outp.pval_vox_p_dpz(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_p_dpz_pervoxel(ivox,ifreq,:))' > d_max_p_dpz)./para.nperm;
        outp.pval_vox_n_dpz(ivox,:,ifreq)   = 1-sum(squeeze(emp.n_n_dpz_pervoxel(ivox,ifreq,:))' > d_max_n_dpz)./para.nperm;
        outp.pval_vox_p_dpz_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_p_context_dpz_pervoxel(ivox,ifreq))' > d_max_ctx_p_dpz)./para.nperm;
        outp.pval_vox_n_dpz_ctx(ivox,ifreq) = 1-sum(squeeze(emp.n_n_context_dpz_pervoxel(ivox,ifreq))' > d_max_ctx_n_dpz)./para.nperm;
      
        outp.pval_vox_p_tvr(ivox,:,ifreq)   = 1-sum(squeeze(emp.taskvsrest_p_pervoxel(ivox,ifreq,:))' > d_max_p_tvr)./para.nperm;
        outp.pval_vox_n_tvr(ivox,:,ifreq)   = 1-sum(squeeze(emp.taskvsrest_n_pervoxel(ivox,ifreq,:))' > d_max_n_tvr)./para.nperm;
      end

    end  
  end
end
