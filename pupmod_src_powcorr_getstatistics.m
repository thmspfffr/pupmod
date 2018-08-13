function [outp, emp] = pupmod_src_powcorr_getstatistics(para)
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
fprintf('Loading data...\n')
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',para.ver));
fprintf('Loading data... Done!\n')
ALPHA = 0.05;
allperms = para.nperm/para.nsubs;

if ~isfield(para,'alpha'); alpha = 0.05; end
if ~isfield(para,'alp'); alp = 0.05; end
if ~isfield(para,'nfreq'); para.nfreq = 13; end
if ~isfield(para,'correction_method'); para.correction_method = 'ranks'; end

if strcmp(para.type,'local') && strcmp(para.correction_method,'single_threshold')
  error('Correction method based on single thresholds is not implemented on a voxel-level! Use ''ranks'' instead.')
end

para.fcsize = size(cleandat,1);
para.nfreq = 13;

% -------------------------------
% Compute number of altered correlations (empirical)
% --------------------------------
emp = pupmod_compute_altered_correlations(cleandat,para);
% --------------------------------

if isfield(para,'stats')
  return 
end
  
clear cleandat s_fc
%% LOAD PERMUTATION DISTRIBUTION
% --------------------------

if strcmp(para.cond,'atx') && strcmp(para.type,'local')
  perm_n_p_atx_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,2,'single');
  perm_n_n_atx_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,2,'single');
elseif strcmp(para.cond,'dpz') && strcmp(para.type,'local')
  perm_n_p_dpz_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,2,'single');
  perm_n_n_dpz_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,2,'single');
elseif strcmp(para.cond,'taskvsrest') && strcmp(para.type,'local')
  perm_taskvsrest_n_p_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,'single');
  perm_taskvsrest_n_n_pervoxel = zeros(para.nperm,para.fcsize,para.nfreq,'single');
end

for iperm = 1 : allperms
  
  fprintf('Load permutation distributions: %d / %d ...\n',iperm,allperms)
  
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,para.nperm,para.ver),'par')
  
  if strcmp(para.cond,'atx')
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
    context_allconn_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_atx_context;
    
  elseif strcmp(para.cond,'dpz')
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
    context_allconn_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_context;
    
  elseif strcmp(para.cond,'taskvsrest')
    % TASK VS REST
    % --------------
    % global effects
    perm_taskvsrest_n_p((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_p;
    perm_taskvsrest_n_n((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_n;
  end
  
  if strcmp(para.cond,'atx') && strcmp(para.type,'local')
    % --------------
    % ATOMOXETINE: local effects
    % --------------
    % rest (indicated by *res1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_res1_pervoxel_n,[2 1 3]);
    % task (indicated by *cnt1* ending)
    perm_n_p_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_p,[2 1 3]);
    perm_n_n_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,2) = permute(par.tperm_cnt1_pervoxel_n,[2 1 3]);
    % context-dependence
    perm_n_p_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_p_pervoxel,[2 1 3]);
    perm_n_n_context_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,1) = permute(par.tperm_context_diff_atx_n_pervoxel,[2 1 3]);
    
  elseif strcmp(para.cond,'dpz') && strcmp(para.type,'local')
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
    
  elseif strcmp(para.cond,'taskvsrest') && strcmp(para.type,'local')
    % local effects
    perm_taskvsrest_n_p_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_p,[2 1 3]);
    perm_taskvsrest_n_n_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = permute(par.tperm_taskvsrest_pervox_n,[2 1 3]);
  end
  
  % double dissociation
  % --------------
  perm_doubledissociation((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_doubledissociation;
end
%%
fprintf('Obtaining corrected p-values ...\n')

if strcmp(para.correction_method, 'ranks')
  % --------------------------
  % MULTIPLE COMPARISONS CORRECTION BASED ON RANKS
  % See Hawellek et al. (2013) Journal of Neuroscience
  % --------------------------
  r = 1 : size(para.nperm,1);
  
  for ifreq = 1 : 13
    
    % GLOBAL EFFECTS
    % ------------------
    if strcmp(para.cond,'atx')
      % ATOMOXETINE
      [~,R_res1_n(:,ifreq)] = sort(perm_n_n_atx(:,ifreq,1),'ascend');
      [~,R_res1_p(:,ifreq)] = sort(perm_n_p_atx(:,ifreq,1),'ascend');
      [~,R_cnt1_n(:,ifreq)] = sort(perm_n_n_atx(:,ifreq,2),'ascend');
      [~,R_cnt1_p(:,ifreq)] = sort(perm_n_p_atx(:,ifreq,2),'ascend');
      
    elseif strcmp(para.cond,'dpz')
      
      [~,R_res2_n(:,ifreq)] = sort(perm_n_n_dpz(:,ifreq,1),'ascend');
      [~,R_res2_p(:,ifreq)] = sort(perm_n_p_dpz(:,ifreq,1),'ascend');
      [~,R_cnt2_n(:,ifreq)] = sort(perm_n_n_dpz(:,ifreq,2),'ascend');
      [~,R_cnt2_p(:,ifreq)] = sort(perm_n_p_dpz(:,ifreq,2),'ascend');
      
    elseif strcmp(para.cond,'taskvsrest')
      % TASK VS REST
      [~,R_tvr_p(:,ifreq)] = sort(perm_taskvsrest_n_p(:,ifreq,1),'ascend');
      [~,R_tvr_n(:,ifreq)] = sort(perm_taskvsrest_n_n(:,ifreq,1),'ascend');
    end
  end
  
  % LOCAL EFFECTS
  % ------------------
  % ATOMOXETINE
  if strcmp(para.cond,'atx') && strcmp(para.type,'local')
    
    D2max_res1_p = squeeze(max(perm_n_p_atx_pervoxel(:,:,:,1),[],2));
    D2max_res1_n = squeeze(max(perm_n_n_atx_pervoxel(:,:,:,1),[],2));
    D2max_cnt1_p = squeeze(max(perm_n_p_atx_pervoxel(:,:,:,2),[],2));
    D2max_cnt1_n = squeeze(max(perm_n_n_atx_pervoxel(:,:,:,2),[],2));
    
    clear perm_n_p_atx_pervoxel perm_n_n_atx_pervoxel 

  % DONEPEZIL
  elseif strcmp(para.cond,'dpz') && strcmp(para.type,'local')
    
    D2max_res2_p = squeeze(max(perm_n_p_dpz_pervoxel(:,:,:,1),[],2));
    D2max_res2_n = squeeze(max(perm_n_n_dpz_pervoxel(:,:,:,1),[],2));
    D2max_cnt2_p = squeeze(max(perm_n_p_dpz_pervoxel(:,:,:,2),[],2));
    D2max_cnt2_n = squeeze(max(perm_n_n_dpz_pervoxel(:,:,:,2),[],2));
    
    clear perm_n_p_dpz_pervoxel perm_n_n_dpz_pervoxel 
   % TASK VS REST
  elseif strcmp(para.cond,'taskvsrest') && strcmp(para.type,'local')
    
    D2max_tvr_p = squeeze(max(perm_taskvsrest_n_p_pervoxel,[],2));
    D2max_tvr_n = squeeze(max(perm_taskvsrest_n_n_pervoxel,[],2));
    
    clear perm_taskvsrest_n_p_pervoxel perm_taskvsrest_n_n_pervoxel 
    
  end
  
  if strcmp(para.cond,'atx')
    Rmax_res1_n = max(R_res1_n,[],2);
    Rmax_res1_p = max(R_res1_p,[],2);
    Rmax_cnt1_p = max(R_cnt1_p,[],2);
    Rmax_cnt1_n = max(R_cnt1_n,[],2);
  elseif strcmp(para.cond,'dpz')
    Rmax_res2_n = max(R_res2_n,[],2);
    Rmax_res2_p = max(R_res2_p,[],2);
    Rmax_cnt2_p = max(R_cnt2_p,[],2);
    Rmax_cnt2_n = max(R_cnt2_n,[],2);
  elseif strcmp(para.cond,'taskvsrest')
    Rmax_tvr_n = max(R_tvr_n,[],2);
    Rmax_tvr_p = max(R_tvr_p,[],2);
  end
  
  if strcmp(para.cond,'atx')
    Dmax_res1_p_corr = zeros(para.nperm,13);
    Dmax_res1_n_corr = zeros(para.nperm,13);
    Dmax_cnt1_p_corr = zeros(para.nperm,13);
    Dmax_cnt1_n_corr = zeros(para.nperm,13);
  elseif strcmp(para.cond,'dpz')
    Dmax_res2_p_corr = zeros(para.nperm,13);
    Dmax_res2_n_corr = zeros(para.nperm,13);
    Dmax_cnt2_p_corr = zeros(para.nperm,13);
    Dmax_cnt2_n_corr = zeros(para.nperm,13);
  elseif strcmp(para.cond,'taskvsrest')
    Dmax_tvr_p_corr = zeros(para.nperm,13);
    Dmax_tvr_n_corr = zeros(para.nperm,13);
  end
  
  for ifreq = 1 : 13
    fprintf('Obtaining corrected p-values: freq%d ...\n',ifreq)
    for irank = 1 : para.nperm
      
      if strcmp(para.cond,'atx')
        Dmax_res1_p_corr(irank,ifreq) = perm_n_p_atx( R_res1_p(:,ifreq) == Rmax_res1_p(irank), ifreq, 1);
        Dmax_res1_n_corr(irank,ifreq) = perm_n_n_atx( R_res1_n(:,ifreq) == Rmax_res1_n(irank), ifreq, 1);
        Dmax_cnt1_p_corr(irank,ifreq) = perm_n_p_atx( R_cnt1_p(:,ifreq) == Rmax_cnt1_p(irank), ifreq, 2);
        Dmax_cnt1_n_corr(irank,ifreq) = perm_n_n_atx( R_cnt1_n(:,ifreq) == Rmax_cnt1_n(irank), ifreq, 2);
      elseif strcmp(para.cond,'dpz')
        Dmax_res2_p_corr(irank,ifreq) = perm_n_p_dpz( R_res2_p(:,ifreq) == Rmax_res2_p(irank), ifreq, 1);
        Dmax_res2_n_corr(irank,ifreq) = perm_n_n_dpz( R_res2_n(:,ifreq) == Rmax_res2_n(irank), ifreq, 1);
        Dmax_cnt2_p_corr(irank,ifreq) = perm_n_p_dpz( R_cnt2_p(:,ifreq) == Rmax_cnt2_p(irank), ifreq, 2);
        Dmax_cnt2_n_corr(irank,ifreq) = perm_n_n_dpz( R_cnt2_n(:,ifreq) == Rmax_cnt2_n(irank), ifreq, 2);
      elseif strcmp(para.cond,'taskvsrest')
        Dmax_tvr_n_corr(irank,ifreq) = perm_taskvsrest_n_p( R_tvr_n(:,ifreq) == Rmax_tvr_n(irank), ifreq);
        Dmax_tvr_p_corr(irank,ifreq) = perm_taskvsrest_n_n( R_tvr_p(:,ifreq) == Rmax_tvr_p(irank), ifreq);
      end
    end
    
    % ---------------------------
    % OBTAIN CORRECTED P-VALUES
    % ---------------------------
    if strcmp(para.cond,'atx')
      outp.p_res1_p(ifreq) = 1-sum(emp.n_p_atx(ifreq,1)>Dmax_res1_p_corr(:,ifreq))/para.nperm;
      outp.p_res1_n(ifreq) = 1-sum(emp.n_n_atx(ifreq,1)>Dmax_res1_n_corr(:,ifreq))/para.nperm;
      outp.p_cnt1_p(ifreq) = 1-sum(emp.n_p_atx(ifreq,2)>Dmax_cnt1_p_corr(:,ifreq))/para.nperm;
      outp.p_cnt1_n(ifreq) = 1-sum(emp.n_n_atx(ifreq,2)>Dmax_cnt1_n_corr(:,ifreq))/para.nperm;
    elseif strcmp(para.cond,'dpz')
      outp.p_res2_p(ifreq) = 1-sum(emp.n_p_dpz(ifreq,1)>Dmax_res2_p_corr(:,ifreq))/para.nperm;
      outp.p_res2_n(ifreq) = 1-sum(emp.n_n_dpz(ifreq,1)>Dmax_res2_n_corr(:,ifreq))/para.nperm;
      outp.p_cnt2_p(ifreq) = 1-sum(emp.n_p_dpz(ifreq,2)>Dmax_cnt2_p_corr(:,ifreq))/para.nperm;
      outp.p_cnt2_n(ifreq) = 1-sum(emp.n_n_dpz(ifreq,2)>Dmax_cnt2_n_corr(:,ifreq))/para.nperm;
    elseif strcmp(para.cond,'taskvsrest')
      outp.p_tvr_p(ifreq) = 1-sum(emp.taskvsrest_p(ifreq)>Dmax_tvr_p_corr(:,ifreq))/para.nperm;
      outp.p_tvr_n(ifreq) = 1-sum(emp.taskvsrest_n(ifreq)>Dmax_tvr_n_corr(:,ifreq))/para.nperm;
    end
  end
  
  % OBTAIN CORRECTED P-VALUES FOR LOCAL TEST (across voxels)
  % ---------------------------
  % significant effects are hard-coded, as the effects are known already.
  % this needs to be changed in a future version of this script, where
  % another permutation test is run after the p-values of the first are
  % known.
  
  if strcmp(para.type,'local')
    if strcmp(para.cond,'atx')
      for irank = 1 : para.nperm
        idx = R_cnt1_p(:,6) == Rmax_cnt1_p(irank);
        D2_adj(irank) = D2max_cnt1_p(idx,6);
      end
      % OBTAIN CORRECTED P-VALUES
      for ivox = 1 : para.fcsize
        outp.p_atx_pervoxel_corrected(ivox) = 1 - sum(emp.n_p_atx_pervoxel(ivox,6,2) > D2_adj) / para.nperm;
      end

    elseif strcmp(para.cond,'dpz')
      for irank = 1 : para.nperm
        idx = R_res2_n(:,7) == Rmax_res2_n(irank);
        D2_adj(irank) = D2max_res2_n(idx,7);
      end
      % OBTAIN CORRECTED P-VALUES
      for ivox = 1 : para.fcsize
        outp.p_dpz_pervoxel_corrected(ivox) = 1 - sum(emp.n_p_dpz_pervoxel(ivox,7,1) > D2_adj) / para.nperm;
      end
    elseif strcmp(para.cond,'taskvsrest')
      error('Does currently not work!')
    end
  end
  
elseif strcmp(para.correction_method,'single_threshold')

  % do not perform rank conversion, but only single threshold test
  % extract maximum fraction of altered corr. across freq first
  if strcmp(para.cond,'atx') && strcmp(para.type,'global')
    % ATOMOXETINE
    % ------------------
    idx_R_res1_p   = max(abs(perm_n_p_atx(:,:,1)),[],2);
    idx_R_res1_n   = max(abs(perm_n_n_atx(:,:,1)),[],2);
    idx_R_cnt1_p   = max(abs(perm_n_p_atx(:,:,2)),[],2);
    idx_R_cnt1_n   = max(abs(perm_n_n_atx(:,:,2)),[],2);
    idx_R_res1_all = max(abs(perm_n_all_atx(:,:,1)),[],2);
    idx_R_cnt1_all = max(abs(perm_n_all_atx(:,:,2)),[],2);
    
  elseif strcmp(para.cond,'dpz') && strcmp(para.type,'global')
    % DONEPEZIL
    % ------------------
    idx_R_res2_p = max(abs(perm_n_p_dpz(:,:,1)),[],2);
    idx_R_res2_n = max(abs(perm_n_n_dpz(:,:,1)),[],2);
    idx_R_cnt2_p = max(abs(perm_n_p_dpz(:,:,2)),[],2);
    idx_R_cnt2_n = max(abs(perm_n_n_dpz(:,:,2)),[],2);
    idx_R_res2_all = max(abs(perm_n_all_dpz(:,:,1)),[],2);
    idx_R_cnt2_all = max(abs(perm_n_all_dpz(:,:,2)),[],2);
    
  elseif strcmp(para.cond,'taskvsrest') && strcmp(para.type,'global')
    % TASK VS REST
    % ------------------
    idx_R_taskvsrest_p        = max(abs(taskvsrest_p_perm),[],2);
    idx_R_taskvsrest_n        = max(abs(taskvsrest_n_perm),[],2);
    
  end
  
  idx_R_doubledissociation  = max(perm_doubledissociation,[],2);
  
  for ifreq = 1 : 13
    if strcmp(para.cond,'atx') && strcmp(para.type,'global')
      % ATOMOXETINE
      % ------------------
      outp.p_res1_p(ifreq) = 1-sum(abs(emp.n_p_atx(ifreq,1))>abs(idx_R_res1_p))/nperm;
      outp.p_res1_n(ifreq) = 1-sum(abs(emp.n_n_atx(ifreq,1))>abs(idx_R_res1_n))/nperm;
      outp.p_cnt1_p(ifreq) = 1-sum(abs(emp.n_p_atx(ifreq,2))>abs(idx_R_cnt1_p))/nperm;
      outp.p_cnt1_n(ifreq) = 1-sum(abs(emp.n_n_atx(ifreq,2))>abs(idx_R_cnt1_n))/nperm;
      outp.p_cnt_atx_all(ifreq) = 1-sum(abs(emp.atx(ifreq,2))>abs(idx_R_cnt1_all))/nperm;
      outp.p_res_atx_all(ifreq) = 1-sum(abs(emp.atx(ifreq,1))>abs(idx_R_res1_all))/nperm;
      
    elseif strcmp(para.cond,'dpz') && strcmp(para.type,'global')
      % DONEPEZIL
      % ------------------
      outp.p_res2_p(ifreq) = 1-sum(abs(emp.n_p_dpz(ifreq,1))>abs(idx_R_res2_p))/nperm;
      outp.p_res2_n(ifreq) = 1-sum(abs(emp.n_n_dpz(ifreq,1))>abs(idx_R_res2_n))/nperm;
      outp.p_cnt2_p(ifreq) = 1-sum(abs(emp.n_p_dpz(ifreq,2))>abs(idx_R_cnt2_p))/nperm;
      outp.p_cnt2_n(ifreq) = 1-sum(abs(emp.n_n_dpz(ifreq,2))>abs(idx_R_cnt2_n))/nperm;
      outp.p_cnt_dpz_all(ifreq) = 1-sum(abs(emp.dpz(ifreq,2))>idx_R_cnt2_all)/nperm;
      outp.p_res_dpz_all(ifreq) = 1-sum(abs(emp.dpz(ifreq,1))>idx_R_res2_all)/nperm;
      
    elseif strcmp(para.cond,'taskvsrest') && strcmp(para.type,'global')
      outp.pval_taskvsrest_p_corr(ifreq) = 1-sum(abs(emp.taskvsrest_p(ifreq))>idx_R_taskvsrest_p(:,ifreq))/nperm;
      outp.pval_taskvsrest_n_corr(ifreq) = 1-sum(abs(emp.taskvsrest_n(ifreq))>idx_R_taskvsrest_n(:,ifreq))/nperm;
    end
  end
end