%% pupmod_all_src_degree
% Computes stats for degree centrality based on a permutation procedure.
% The pharma labels are exchanged randomly *nperm* times and degree is
% computed for each permutation. For each permutation, the maximum degree
% across space is noted, resulting in a maximum degree distribution Dmax.
% From Dmax, p-values can be obtained.

clear

v = 12;

outdir = '~/pupmod/proc/conn/';

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));


%%

s_res = squeeze(cleandat(:,:,:,:,1,:));
s_cnt = squeeze(cleandat(:,:,:,:,2,:));

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')

nperm = 10000; alp = 0.05;

par.subs = 50;
par.allperms = nperm/par.subs;

if ~exist(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v),'all_idx1');
else
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v));
end

dat_cnt1 = squeeze(cleandat(:,:,:,[1 2],2,:));
dat_res1 = squeeze(cleandat(:,:,:,[1 2],1,:));
dat_cnt2 = squeeze(cleandat(:,:,:,[1 3],2,:));
dat_res2 = squeeze(cleandat(:,:,:,[1 3],1,:));

for iperm = 1 : par.allperms
  
  if ~exist(sprintf(['~/pupmod/proc/' 'pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' '~/pupmod/proc/' sprintf('pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    
    
    for icond = 1 : 2
      for ifoi = 1:13
        
        if ~exist(sprintf([outdir 'pupmod_all_src_degree_permtest_c%d_f%d_v%d_processing.txt'],icond,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_all_src_permtest_degree_c%d_f%d_v%d_processing.txt',icond,ifoi,v)]);
        else
          continue
        end
        
        j = 1:size(cleandat,1);
        
        for il = 1 : size(cleandat,1)
          
          jj = j(j~=il);
          
          for jl = 1 : size(cleandat,1)
            
            jjj = j(j~=jl);
            
            fc_tmp = cleandat(:,:,:,1:3,icond,ifoi);
            
            fprintf('Computing location %d %d ...\n',il,jl);
            
            x = squeeze(fc_tmp(il,jl,:,:));
            
            x_ref1 = squeeze(nanmean(fc_tmp(il,jj,:,:),4));
            x_ref2 = squeeze(nanmean(fc_tmp(jjj,jl,:,:),4));
            
            m1 = mean(x_ref1); s1 = std(x_ref1);
            m2 = mean(x_ref2); s2 = std(x_ref2);
            
            z11 = (x(:,1)' - m1)./s1;
            z12 = (x(:,1)' - m2)./s2;
            
            z21 = (x(:,2)' - m1)./s1;
            z22 = (x(:,2)' - m2)./s2;
            
            z31 = (x(:,3)' - m1)./s1;
            z32 = (x(:,3)' - m2)./s2;
            
            [~,p11] = ttest(z11);
            [~,p12] = ttest(z12);
            
            [~,p21] = ttest(z21);
            [~,p22] = ttest(z22);
            
            [~,p31] = ttest(z31);
            [~,p32] = ttest(z32);
            
            th11(jl,:) = p11 < 0.05/2;
            th12(jl,:) = p12 < 0.05/2;
            th21(jl,:) = p21 < 0.05/2;
            th22(jl,:) = p22 < 0.05/2;
            th31(jl,:) = p31 < 0.05/2;
            th32(jl,:) = p32 < 0.05/2;
            
          end
          
          th11 = th11(jj,:);
          th12 = th12(jj,:);
          th21 = th21(jj,:);
          th22 = th22(jj,:);
          th31 = th31(jj,:);
          th32 = th32(jj,:);
          
          % any connection significant?
          th1(il,jj,:) = (th11 + th12) > 0; clear th11 th12
          th2(il,jj,:) = (th21 + th22) > 0; clear th21 th22
          th3(il,jj,:) = (th31 + th32) > 0; clear th31 th32
          
        end
        
        k(:,1) = sum(th1)/size(cleandat,1);
        k(:,2) = sum(th2)/size(cleandat,1);
        para.type = 'global';
        para.nperm = 50000;
        para.nsubs = 100;
        para.ver = 12;
        para.correction_method = 'ranks';
        
        para.cond = 'taskvsrest';
        
        p = pupmod_src_powcorr_getstatistics(para);
        k(:,3) = sum(th3)/size(cleandat,1);
        
        save(sprintf([outdir 'pupmod_all_src_degree_permtest_c%d_f%d_v%d.mat'],icond,ifoi,v),'k');
        
      end
    end
    
    error('!')
    
    % Add more metrics
    
   