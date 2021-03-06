%% pupmod_all_src_degree_permtest
% Computes stats for degree centrality based on a permutation procedure.
% The pharma labels are exchanged randomly *nperm* times and degree is
% computed for each permutation. For each permutation, the maximum degree
% across space is noted, resulting in a maximum degree distribution Dmax.
% From Dmax, p-values can be obtained.

clear
% -----------------
% VERSION 12 - see pupmod_src_powcorr.m
% -----------------
v = 23;
nperm = 10000; 
par.subs = 250;
par.allperms = nperm/par.subs;
% -----------------
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

%%

mask    = logical(triu(ones(400,400),1));
para = [];
para.alpha = 0.01;
para.nfreq = 25;
para.absolute = 0;
para.relative_degree = 1;
para.clustering = 1;
para.transitivity = 1;
para.charpath = 1;

cleandat = pupmod_loadpowcorr(v,SUBJLIST,1);
fcsize = size(cleandat,1);

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')


if ~exist(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v),'all_idx1');
else
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v));
end

% dat_atx     = squeeze(cleandat(:,:,:,[1 2],:,:));
% dat_dpz     = squeeze(cleandat(:,:,:,[1 3],:,:));
% taskvsrest  = squeeze(cleandat(:,:,:,1,:,:));

fcsize = size(cleandat,1);

% clear cleandat

for iperm = 1 : par.allperms
%   
  if ~exist(sprintf(['~/pupmod/proc/' 'pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' '~/pupmod/proc/' sprintf('pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
%   
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    
    fprintf('Perm #%d ...\n',kperm);
    
    % -----------------------------
    % RESHUFFLE THE LABELS (within subjects)
    % -----------------------------
    % Atomoxetine
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat_atx(:,:,i,1,:,:) = cleandat(:,:,i,idx1(i),:,:);
      permdat_atx(:,:,i,2,:,:) = cleandat(:,:,i,idx2(i),:,:);
    end
    
    % Donepezil
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    idx1(idx1==2) = 3;
    idx2(idx2==2) = 3;
    
    for i = 1 : length(idx1)
      permdat_dpz(:,:,i,1,:,:) = cleandat(:,:,i,idx1(i),:,:);
      permdat_dpz(:,:,i,2,:,:) = cleandat(:,:,i,idx2(i),:,:);
    end
    
    % Task vs rest
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      taskvsrest_perm(:,:,i,1,:) = cleandat(:,:,i,1,idx1(i),:);
      taskvsrest_perm(:,:,i,2,:) = cleandat(:,:,i,1,idx2(i),:);
    end

    for icond = 1 : 2
        
%       para = [];
      para.alpha = 0.01;
      para.nfreq = 25;
      
      graph_atx = tp_degree(permdat_atx(:,:,:,:,icond,:),para);
      graph_dpz = tp_degree(permdat_atx(:,:,:,:,icond,:),para);

      % total degree
      outp.perm_k_atx(:,:,icond,kperm) = graph_atx.tot_degree;
      outp.perm_k_dpz(:,:,icond,kperm) = graph_dpz.tot_degree;
      
      % degree per voxel
      outp.perm_k_atx_pervoxel(:,:,:,icond,kperm) = graph_atx.node_degree;
      outp.perm_k_dpz_pervoxel(:,:,:,icond,kperm) = graph_dpz.node_degree;
            
      % clustering coeff
      outp.perm_clust_atx(:,:,icond,kperm) = squeeze(mean(graph_atx.clust));
      outp.perm_clust_dpz(:,:,icond,kperm) = squeeze(mean(graph_dpz.clust));
      
      % clustering coeff per voxel
      outp.perm_clust_atx_pervoxel(:,:,:,icond,kperm) = graph_atx.clust;
      outp.perm_clust_dpz_pervoxel(:,:,:,icond,kperm) = graph_dpz.clust;
      
      % transitivity
      outp.perm_trans_atx(:,:,icond,kperm) = graph_atx.transitivity;
      outp.perm_trans_dpz(:,:,icond,kperm) = graph_dpz.transitivity;
      
      % char path
      outp.perm_charpath_atx(:,:,icond,kperm) = graph_atx.charpath;
      outp.perm_charpath_dpz(:,:,icond,kperm) = graph_dpz.charpath;
      
      % compute both for task vs rest
%       deg_tvr = tp_degree(taskvsrest_perm,para);
%       outp.perm_k_tvr(:,:,icond,kperm) = nansum(reshape(deg_tvr,[fcsize^2 para.nfreq 2]))/fcsize^2;
%       outp.perm_k_tvr_pervoxel(:,:,icond,kperm) = squeeze(nansum((deg_tvr)))/fcsize;
    
%       deg_atx(:,:,:,:,icond) = tp_degree(permdat_atx(:,:,:,:,icond,:),para);
%       deg_dpz(:,:,:,:,icond) = tp_degree(permdat_dpz(:,:,:,:,icond,:),para);

    end
    
    % compute global degree 
   
  end
    
  save(sprintf(['~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'outp');

end

error('!')

%% CONCATENATE PERMUTATONS AND COMPUTE STATS

% for iperm = 1 : par.allperms
%   
%   load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
%   
%   perm_k_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = outp.perm_k_atx;
%   perm_k_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = outp.perm_k_dpz;
%   perm_k_atx_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,:) = outp.perm_k_atx_pervoxel;
%   perm_k_dpz_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,:) = outp.perm_k_dpz_pervoxel;
%   perm_k_tvr((iperm-1)*par.subs+1:(iperm)*par.subs,:,:) = outp.perm_k_tvr;
%   perm_k_tvr_pervoxel((iperm-1)*par.subs+1:(iperm)*par.subs,:,:,:) = outp.perm_k_tvr_pervoxel;
%   
% end
