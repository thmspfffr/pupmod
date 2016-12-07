%% pupmod_src_powcorr_permtest

clear all

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/


%%
clear s s1 s2 fc_mean
v = 10;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        s_cnt(:,:,isubj,m,ifoi,iblock) =  powcorr;
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        s_res(:,:,isubj,m,ifoi,iblock) =  powcorr;
        
      end
    end
  end
end

% save('~/pupmod/proc/pupmod_src_fcd.mat','h')

s_res = squeeze(nanmean(s_res(:,:,SUBJLIST,:,:,:),6));
s_cnt = squeeze(nanmean(s_cnt(:,:,SUBJLIST,:,:,:),6));

%% PERMUTATION TEST
% nperm = 1000;
%
% rng('shuffle');
%
% dat   = s_cnt(:,:,:,[1 3],:);
%
% for iperm = 1 : nperm
%
%
%   % within subjects permutation test
%
%   disp(sprintf('Perm #%d',iperm));
%
%  	idx1 = randi(2,[size(s_cnt,3),1]);
% 	idx2 = 3-idx1;
%
%   for i = 1 : length(idx1)
%
%     x(:,:,i,1,:) = dat(:,:,i,idx1(i),:);
%     x(:,:,i,2,:) = dat(:,:,i,idx2(i),:);
%
%   end
%
%   for ifoi = 1 : 13
%
%     [~,~,~,stat] = ttest(x(:,:,:,2,ifoi),x(:,:,:,1,ifoi),'dim',3);
%     tmax(iperm,ifoi) = max(stat.tstat(:));
%
%   end
%
% end

%%
v = 10;
clear permdat_cnt1
clear permdat_cnt2
clear permdat_res1
clear permdat_res2

rng('shuffle')

nperm = 10000; alp = 0.2;

par.subs = 200;
par.allperms = nperm/par.subs;

if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v),'all_idx1');
else
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v));
end

dat_cnt1 = s_cnt(:,:,:,[1 2],:);
dat_res1 = s_res(:,:,:,[1 2],:);
dat_cnt2 = s_cnt(:,:,:,[1 3],:);
dat_res2 = s_res(:,:,:,[1 3],:);

for iperm = 1 : par.allperms
  
  if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' outdir sprintf('pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    
    % within subjects permutation test
    
    disp(sprintf('Perm #%d',kperm));
    
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
      permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);
      
      permdat_res1(:,:,i,1,:) = dat_res1(:,:,i,idx1(i),:);
      permdat_res1(:,:,i,2,:) = dat_res1(:,:,i,idx2(i),:);
      
    end
    
    idx1 = randi(2,[size(dat_cnt1,3),1]);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt2(:,:,i,1,:) = dat_cnt2(:,:,i,idx1(i),:);
      permdat_cnt2(:,:,i,2,:) = dat_cnt2(:,:,i,idx2(i),:);
      
      permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
      permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:);
      
    end
    
    for ifoi = 1 : 13
      
      [t_cnt1,~,~,s] = ttest(permdat_cnt1(:,:,:,2,ifoi),permdat_cnt1(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_cnt1 = t_cnt1.*sign(s.tstat); clear s
      [t_res1,~,~,s] = ttest(permdat_res1(:,:,:,2,ifoi),permdat_res1(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_res1 = t_res1.*sign(s.tstat); clear s
      [t_cnt2,~,~,s] = ttest(permdat_cnt2(:,:,:,2,ifoi),permdat_cnt2(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_cnt2 = t_cnt2.*sign(s.tstat); clear s
      [t_res2,~,~,s] = ttest(permdat_res2(:,:,:,2,ifoi),permdat_res2(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_res2 = t_res2.*sign(s.tstat); clear s
      
      [t_all1,~,~,s] = ttest(permdat_res1(:,:,:,2,ifoi)-permdat_res1(:,:,:,1,ifoi),permdat_cnt1(:,:,:,2,ifoi)-permdat_cnt1(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_all1 = t_all1.*sign(s.tstat); clear s
      [t_all2,~,~,s] = ttest(permdat_res2(:,:,:,2,ifoi)-permdat_res2(:,:,:,1,ifoi),permdat_cnt2(:,:,:,2,ifoi)-permdat_cnt2(:,:,:,1,ifoi),'dim',3,'alpha',alp);
      t_all2 = t_all2.*sign(s.tstat); clear s
      
      par.tperm_cnt1_n(kperm,ifoi)=sum(sum(triu(t_cnt1<0,1)));
      par.tperm_cnt1_p(kperm,ifoi)=sum(sum(triu(t_cnt1>0,1)));
      
      par.tperm_res1_n(kperm,ifoi)=sum(sum(triu(t_res1<0,1)));
      par.tperm_res1_p(kperm,ifoi)=sum(sum(triu(t_res1>0,1)));
      
      par.tperm_cnt2_n(kperm,ifoi)=sum(sum(triu(t_cnt2<0,1)));
      par.tperm_cnt2_p(kperm,ifoi)=sum(sum(triu(t_cnt2>0,1)));
      
      par.tperm_res2_n(kperm,ifoi)=sum(sum(triu(t_res2<0,1)));
      par.tperm_res2_p(kperm,ifoi)=sum(sum(triu(t_res2>0,1)));
      
      par.tperm_all1_n(kperm,ifoi)=sum(sum(triu(t_all1<0,1)));
      par.tperm_all1_p(kperm,ifoi)=sum(sum(triu(t_all1>0,1)));
      
      par.tperm_all2_n(kperm,ifoi)=sum(sum(triu(t_all2<0,1)));
      par.tperm_all2_p(kperm,ifoi)=sum(sum(triu(t_all2>0,1)));
      
    end
  end
  
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
  
  try
    pause(randi(3))
    load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  catch me
    save(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
  end
  
end
%%
