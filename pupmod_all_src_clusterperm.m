%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 1;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

ord = pconn_randomization;

    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));


for ifoi = 1:13
%     ifoi
    
    s_fc(:,:,:,:,1) = cleandat(:,:,:,:,1,ifoi);
    s_fc(:,:,:,:,2) = cleandat(:,:,:,:,2,ifoi);
    
%     for isubj = SUBJLIST
%         disp(isubj)
%         for m = 1 : 3
%             
%             im = find(ord(isubj,:)==m);
%             
%             for iblock = 1 : 2
%                 clear tmp
%                 load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%                 
%                 p1(:,:,iblock) = single(powcorr); clear powcorr
%                 
%                 load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%                 
%                 p2(:,:,iblock) = single(powcorr); clear powcorr
%                                 
%             end
%             
%             s_fc(:,:,isubj,m,1) = nanmean(p1,3);
%             s_fc(:,:,isubj,m,2) = nanmean(p2,3);
%             
%             clear p1 p2
%                        
%         end
%     end
%     
%     s_fc = s_fc(:,:,SUBJLIST,:,:,:);

    %
    
    [h,~,~,s]=ttest(s_fc(:,:,:,2,1),s_fc(:,:,:,1,1),'dim',3);
    n_p_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
        
    [h,~,~,s]=ttest(s_fc(:,:,:,2,2),s_fc(:,:,:,1,2),'dim',3);
    n_p_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  
    [h,~,~,s]=ttest(s_fc(:,:,:,3,1),s_fc(:,:,:,1,1),'dim',3);
    n_p_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    
    [h,~,~,s]=ttest(s_fc(:,:,:,3,2),s_fc(:,:,:,1,2),'dim',3);
    n_p_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    [h,~,~,s]=ttest(s_fc(:,:,:,2,1)-s_fc(:,:,:,1,1),s_fc(:,:,:,2,2)-s_fc(:,:,:,1,2),'dim',3);
    n_p_dd_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_dd_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    context_allconn_emp_atx_test(ifoi) = nansum(nansum((abs(h))))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    [h,~,~,s]=ttest(s_fc(:,:,:,3,1)-s_fc(:,:,:,1,1),s_fc(:,:,:,3,2)-s_fc(:,:,:,1,2),'dim',3);
    n_p_dd_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    n_n_dd_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    %   context_allconn_emp_dpz(ifoi) = nansum(nansum((abs(h))))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    context_allconn_emp_dpz_test(ifoi) = nansum(nansum((abs(h))))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    %   NUMBER OF ALTERED CONNECTIONS
    % --------------------------
    h=ttest(s_fc(:,:,:,2,1),s_fc(:,:,:,1,1),'dim',3);
    atx(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    h=ttest(s_fc(:,:,:,2,2),s_fc(:,:,:,1,2),'dim',3);
    atx(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    h=ttest(s_fc(:,:,:,3,1),s_fc(:,:,:,1,1),'dim',3);
    dpz(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    h=ttest(s_fc(:,:,:,3,2),s_fc(:,:,:,1,2),'dim',3);
    dpz(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
    
    context_allconn_emp_atx(ifoi) = atx(ifoi,2)-atx(ifoi,1);
    context_allconn_emp_dpz(ifoi) = dpz(ifoi,2)-dpz(ifoi,1);
    
    
end

clear tperm_cnt1_n tperm_cnt1_p tperm_cnt2_n tperm_cnt2_p
clear tperm_res1_n tperm_res1_p tperm_res2_n tperm_res2_p
clear perm_n_p_atx perm_n_n_atx perm_n_n_dpz perm_n_p_dpz

nperm = 20000;
v = 1;
par.subs = 100;
par.allperms = nperm/par.subs;
alpha = 0.05;
alp = 0.05;

for iperm = 1 : par.allperms
    
    load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
    
    perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_p;
    perm_n_p_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt1_p;
    
    perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res1_n;
    perm_n_n_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt1_n;
    
    perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_p;
    perm_n_p_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt2_p;
    
    perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_res2_n;
    perm_n_n_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_cnt2_n;
    
    perm_n_all_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_atx_during_rest;
    perm_n_all_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_atx_during_task;
    
    perm_n_all_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_during_rest;
    perm_n_all_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,2)=par.tperm_dpz_during_task;
    
    try
        context_allconn_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.t;
    catch me
        context_allconn_atx((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_atx_context;
    end
    
    context_allconn_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_context;
    
    context_allconn_atx_test((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_atx_context_test;
    context_allconn_dpz_test((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_context_test;
    
end


for ifoi = 1 : 13
    
    % COMPUTE EFFECT OF DRUGS ON TASK
    tp_cnt1_n(ifoi) = 1-[sum(n_n_atx(ifoi,2)>perm_n_n_atx(:,ifoi,2))]/nperm;
    tp_cnt1_p(ifoi) = 1-[sum(n_p_atx(ifoi,2)>perm_n_p_atx(:,ifoi,2))]/nperm;
    %
    tp_cnt2_n(ifoi) = 1-[sum(n_n_dpz(ifoi,2)>perm_n_n_dpz(:,ifoi,2))]/nperm;
    tp_cnt2_p(ifoi) = 1-[sum(n_p_dpz(ifoi,2)>perm_n_p_dpz(:,ifoi,2))]/nperm;
    
    % COMPUTE EFFECT OF DRUGS ON REST
    tp_res1_n(ifoi) = 1-[sum(n_n_atx(ifoi,1)>perm_n_n_atx(:,ifoi,1))]/nperm;
    tp_res1_p(ifoi) = 1-[sum(n_p_atx(ifoi,1)>perm_n_p_atx(:,ifoi,1))]/nperm;
    
    tp_res2_n(ifoi) = 1-[sum(n_n_dpz(ifoi,1)>perm_n_n_dpz(:,ifoi,1))]/nperm;
    tp_res2_p(ifoi) = 1-[sum(n_p_dpz(ifoi,1)>perm_n_p_dpz(:,ifoi,1))]/nperm;
    
    % COMPUTE EFFECT ON ALL CONNECTIONS
    tp_res1_all(ifoi) = 1-[sum(atx(ifoi,1)>perm_n_all_atx(:,ifoi,1))]/nperm;
    tp_res2_all(ifoi) = 1-[sum(dpz(ifoi,1)>perm_n_all_dpz(:,ifoi,1))]/nperm;
    
    tp_cnt1_all(ifoi) = 1-[sum(atx(ifoi,2)>perm_n_all_atx(:,ifoi,2))]/nperm;
    tp_cnt2_all(ifoi) = 1-[sum(dpz(ifoi,2)>perm_n_all_dpz(:,ifoi,2))]/nperm;
    
    %   % COMPUTE CHANGE OF CHANGE (no need to test two-sided)
    p_d1_p(ifoi) = 1-sum(n_p_atx(ifoi,1)-n_p_atx(ifoi,2)<(perm_n_p_atx(:,ifoi,1)-perm_n_p_atx(:,ifoi,2)))/nperm;
    p_d1_n(ifoi) = 1-sum(abs(n_n_atx(ifoi,1)-n_n_atx(ifoi,2))>abs((perm_n_n_atx(:,ifoi,1)-perm_n_n_atx(:,ifoi,2))))/nperm;
    %
    p_d2_p(ifoi) = 1-sum(abs(n_p_dpz(ifoi,1)-n_p_dpz(ifoi,2))>abs((perm_n_p_dpz(:,ifoi,1)-perm_n_p_dpz(:,ifoi,2))))/nperm;
    p_d2_n(ifoi) = 1-sum(n_n_dpz(ifoi,1)-n_n_dpz(ifoi,2)>(perm_n_n_dpz(:,ifoi,1)-perm_n_n_dpz(:,ifoi,2)))/nperm;
    
    % this is the version where the difference in the number of
    % significantly altered connections is tested between rest and task
    pval_context_atx_all(ifoi) = 1-[sum(abs(context_allconn_emp_atx(ifoi))>context_allconn_atx(:,ifoi))]/nperm;
    pval_context_dpz_all(ifoi) = 1-[sum(context_allconn_emp_dpz(ifoi)<context_allconn_dpz(:,ifoi))]/nperm;
    
    % this is the version where the number of context-dependent connections
    % is counted and the p-value is derived
    pval_context_atx_all_test(ifoi) = 1-[sum(abs(context_allconn_emp_atx_test(ifoi))>context_allconn_atx_test(:,ifoi))]/nperm;
    pval_context_dpz_all_test(ifoi) = 1-[sum(context_allconn_emp_dpz_test(ifoi)<context_allconn_dpz_test(:,ifoi))]/nperm;
    
end

n_clust = max(bwlabel(tp_cnt1_p<0.05));

cnt = 0;
for iclust = 1 : n_clust  
  if length(find(bwlabel(tp_cnt1_p<0.05)==iclust)) < 2    
    continue   
  else
    cnt = cnt + 1;
    emp_clust_sum(cnt) = sum(tinv(tp_cnt1_p(find(bwlabel(tp_cnt1_p<0.05)==iclust)),28));
    
  end
end
   
clear s_fc
%%
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

s_res = squeeze(cleandat(:,:,:,:,1,:));
s_cnt = squeeze(cleandat(:,:,:,:,2,:));

dat_cnt1 = s_cnt(:,:,:,[1 2],:);  
dat_res1 = s_res(:,:,:,[1 2],:);  
dat_cnt2 = s_cnt(:,:,:,[1 3],:); clear s_cnt
dat_res2 = s_res(:,:,:,[1 3],:); clear s_res

for iperm = 1 : 500
 
    disp(sprintf('Perm #%d',iperm));
    
    idx1 = randi(2,[size(SUBJLIST,2),1]);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
      permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);
      
%       permdat_res1(:,:,i,1,:) = dat_res1(:,:,i,idx1(i),:);
%       permdat_res1(:,:,i,2,:) = dat_res1(:,:,i,idx2(i),:);
      
    end
    
%     for i = 1 : length(idx1)
%       
%       permdat_cnt2(:,:,i,1,:) = dat_cnt2(:,:,i,idx1(i),:);
%       permdat_cnt2(:,:,i,2,:) = dat_cnt2(:,:,i,idx2(i),:);
%       
%       permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
%       permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:);
%       
%     end
    
    
      
      % compute ttest during task and atomoxetine
      [t_cnt1,~,~,s] = ttest(permdat_cnt1(:,:,:,2,:),permdat_cnt1(:,:,:,1,:),'dim',3,'alpha',alp);
     t_cnt1 = squeeze(t_cnt1);
     s.tstat = squeeze(s.tstat);
      t_cnt1 = t_cnt1.*sign(s.tstat); 

        
      % compute ttest during rest and atomoxetine
%       [t_res1,~,~,s] = ttest(permdat_res1(:,:,:,2,ifoi),permdat_res1(:,:,:,1,ifoi),'dim',3,'alpha',alp);
%       t_res1 = t_res1.*sign(s.tstat); clear s
      tperm_cnt1_p=squeeze(nansum(nansum(t_cnt1>0)))./8010;
      
      for ifoi = 1 : 13
        tp_cnt1_p(ifoi) = 1-[sum(tperm_cnt1_p(ifoi)>perm_n_p_atx(:,ifoi,2))]/nperm;
      

    end
    clear s
    
    n_clust = max(bwlabel(tp_cnt1_p<0.05));
    cnt = 0;
    if n_clust == 0 
      clust_sum = tinv(min(tp_cnt1_p),27);
    elseif n_clust > 0 
      for iclust = 1 : n_clust  
        if length(find(bwlabel(tp_cnt1_p<0.05)==iclust)) < 2  
          clust_sum = tinv(tp_cnt1_p(find(bwlabel(tp_cnt1_p<0.05)==iclust)),27);
          continue   
        else
          cnt = cnt + 1;
          clust_sum(cnt) = sum(tinv(tp_cnt1_p(find(bwlabel(tp_cnt1_p<0.05)==iclust)),27));    
        end
      end
    else
      warning('warning')
    end
    
  

    max_clust(iperm) = max(clust_sum)
    if max(clust_sum) == 0
      error('!')
    end
      
         clear clust_sum
 
    
end



