%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 14;

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
    
    % DOUBLE DISSOCIATION ACROSS FREUQUENCIES -----------------------------
    doubledissociation_emp = context_allconn_emp_atx-context_allconn_emp_dpz;
    
    % TASK VS REST ------------------------------------------------------
    [t_tvsr,~,~,s] = ttest(s_fc(:,:,:,1,2),s_fc(:,:,:,1,1),'dim',3);
    t_tvsr = t_tvsr.*sign(s.tstat); clear s   
      
    taskvsrest_p(ifoi) = nansum(nansum(t_tvsr>0))./8010;
    taskvsrest_n(ifoi) = nansum(nansum(t_tvsr<0))./8010;
      
end

error('!')

%% PLOT NUMBER OF ALTERED CONNECTIONS (IRRESPECTIVE OF SIGN)

figure;

subplot(1,2,1); hold on

plot(atx(:,1),'color',[1 0.7 0.4],'linewidth',2)
plot(atx(:,2),'color',[1 0.1 0.1],'linewidth',2)

axis([0 14 -0.02 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

xlabel('Frequency [Hz]'); ylabel('Fraction of alt. corr. [%]')
axis square
subplot(1,2,2); hold on

plot(dpz(:,1),'color',[0.4 0.7 1],'linewidth',2)
plot(dpz(:,2),'color',[0.1 0.1 1],'linewidth',2)

axis([0 14 -0.02 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

xlabel('Frequency [Hz]'); ylabel('Fraction of alt. corr. [%]')

axis square
%% PLOT

figure;

subplot(3,2,1); hold on

plot(n_p_atx(:,1),'r')
plot(n_n_atx(:,1),'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,3); hold on

plot(n_p_atx(:,2),'r')
plot(n_n_atx(:,2),'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,2); hold on

plot(n_p_dpz(:,1),'r')
plot(n_n_dpz(:,1),'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,4); hold on

plot(n_p_dpz(:,2),'r')
plot(n_n_dpz(:,2),'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,5); hold on

plot(n_p_dd_atx,'r')
plot(n_n_dd_atx(:,1),'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,6); hold on

plot(n_p_dd_dpz,'r')
plot(n_n_dd_dpz,'b')

axis([0 14 0 0.5])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

%%
clear tperm_cnt1_n tperm_cnt1_p tperm_cnt2_n tperm_cnt2_p
clear tperm_res1_n tperm_res1_p tperm_res2_n tperm_res2_p
clear perm_n_p_atx perm_n_n_atx perm_n_n_dpz perm_n_p_dpz

nperm = 10000;
v = 14;
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
    
    perm_doubledissociation((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_doubledissociation;
    
    context_allconn_dpz((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_context;
    
    context_allconn_atx_test((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_atx_context_test;
    context_allconn_dpz_test((iperm-1)*par.subs+1:(iperm)*par.subs,:,1)=par.tperm_dpz_context_test;
    
    taskvsrest_p_perm((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_p;
    taskvsrest_n_perm((iperm-1)*par.subs+1:(iperm)*par.subs,:) = par.tperm_taskvsrest_n;

end



%%

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
    
    % COMPUTe EFFECT ON ALL CONNECTIONS
    tp_res1_all(ifoi) = 1-[sum(atx(ifoi,1)>perm_n_all_atx(:,ifoi,1))]/nperm;
    tp_res2_all(ifoi) = 1-[sum(dpz(ifoi,1)>perm_n_all_dpz(:,ifoi,1))]/nperm;
    
    tp_cnt1_all(ifoi) = 1-[sum(atx(ifoi,2)>perm_n_all_atx(:,ifoi,2))]/nperm;
    tp_cnt2_all(ifoi) = 1-[sum(dpz(ifoi,2)>perm_n_all_dpz(:,ifoi,2))]/nperm;
    
    %   % COMPUTE CHANGE OF CHANGE (no need to test two-sided)
    p_d1_p(ifoi) = 1-sum(abs(n_p_atx(ifoi,1)-n_p_atx(ifoi,2))>abs((perm_n_p_atx(:,ifoi,1)-perm_n_p_atx(:,ifoi,2))))/nperm;
    p_d1_n(ifoi) = 1-sum(abs(n_n_atx(ifoi,1)-n_n_atx(ifoi,2))>abs((perm_n_n_atx(:,ifoi,1)-perm_n_n_atx(:,ifoi,2))))/nperm;
    %
    p_d2_p(ifoi) = 1-sum(abs(n_p_dpz(ifoi,1)-n_p_dpz(ifoi,2))>abs((perm_n_p_dpz(:,ifoi,1)-perm_n_p_dpz(:,ifoi,2))))/nperm;
    p_d2_n(ifoi) = 1-sum(abs(n_n_dpz(ifoi,1)-n_n_dpz(ifoi,2))>abs((perm_n_n_dpz(:,ifoi,1)-perm_n_n_dpz(:,ifoi,2))))/nperm;
    
    % this is the version where the difference in the number of
    % significantly altered connections is tested between rest and task
    pval_context_atx_all(ifoi) = 1-[sum(abs(context_allconn_emp_atx(ifoi))>abs(context_allconn_atx(:,ifoi)))]/nperm;
    pval_context_dpz_all(ifoi) = 1-[sum(abs(context_allconn_emp_dpz(ifoi))>abs(context_allconn_dpz(:,ifoi)))]/nperm;
    
    % this is the version where the number of context-dependent connections
    % is counted and the p-value is derived
    pval_context_atx_all_test(ifoi) = 1-[sum(abs(context_allconn_emp_atx_test(ifoi))>abs(context_allconn_atx_test(:,ifoi)))]/nperm;
    pval_context_dpz_all_test(ifoi) = 1-[sum(abs(context_allconn_emp_dpz_test(ifoi))>abs(context_allconn_dpz_test(:,ifoi)))]/nperm;

    pval_doubledissociation_all(ifoi) = 1-[sum(abs(doubledissociation_emp(ifoi))>abs(perm_doubledissociation(:,ifoi)))]/nperm;
    
    pval_taskvsrest_p(ifoi) = 1-[sum(abs(taskvsrest_p(ifoi))>abs(taskvsrest_p_perm(:,ifoi)))]/nperm;
    pval_taskvsrest_n(ifoi) = 1-[sum(abs(taskvsrest_n(ifoi))>abs(taskvsrest_n_perm(:,ifoi)))]/nperm;
    
end


%% PLOT NUMBER OF ALTERED CONNECTIONS (IRRESPECTIVE OF SIGN)

figure;

subplot(2,2,1); hold on

plot(atx(:,1).*100,'color',[1 0.7 0.4],'linewidth',2)
plot(atx(:,2).*100,'color',[1 0.1 0.1],'linewidth',2)

plot(prctile(perm_n_all_atx(:,:,1).*100,95),'linewidth',1,'color',[1 0.5 0.2],'linestyle',':')
plot(prctile(perm_n_all_atx(:,:,2).*100,95),'linewidth',1,'color',[1 0.1 0.1],'linestyle',':')

axis([0 14 -2 30])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

xlabel('Frequency [Hz]'); ylabel('Fraction of alt. corr. [%]')
axis square
subplot(2,2,2); hold on

plot(dpz(:,1).*100,'color',[0.4 0.7 1],'linewidth',2)
plot(dpz(:,2).*100,'color',[0.1 0.1 1],'linewidth',2)

plot(prctile(perm_n_all_dpz(:,:,1).*100,95),'linewidth',1,'color',[0.4 0.7 1],'linestyle',':')
plot(prctile(perm_n_all_dpz(:,:,2).*100,95),'linewidth',1,'color',[0.1 0.1 1],'linestyle',':')


axis([0 14 -2 30])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

xlabel('Frequency [Hz]'); ylabel('Fraction of alt. corr. [%]')

axis square
print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_alteredcorrelations.pdf'))


%% PLOT CHANGES

figure;

subplot(3,2,1); hold on
plot(n_p_atx(:,1),'r-','linewidth',3)
plot(n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3],'yticklabel',[0 10 20 30])
ylabel('Altered corr. [%]')
title('Rest')
axis([1 14 -0.05 0.55])
plot(find(tp_res1_p<0.025),n_p_atx(find(tp_res1_p<0.025),1),'k.','markersize',30)
plot(find(tp_res1_n<0.025),n_n_atx(find(tp_res1_n<0.025),1),'k.','markersize',30)

plot(prctile(perm_n_n_atx(:,:,1),95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(perm_n_p_atx(:,:,1),95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,2); hold on
plot(n_p_dpz(:,1),'r-','linewidth',3)
plot(n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([1 14 -0.05 0.55])
plot(find(tp_res2_p<0.025),n_p_dpz(find(tp_res2_p<0.025),1),'k.','markersize',30)
plot(find(tp_res2_n<0.025),n_n_dpz(find(tp_res2_n<0.025),1),'k.','markersize',30)

plot(prctile(perm_n_n_dpz(:,:,1),95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(perm_n_p_dpz(:,:,1),95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,3); hold on
plot(n_p_atx(:,2),'r-','linewidth',3)
plot(n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
ylabel('Altered corr. [%]')
title('Task')
axis([1 14 -0.05 0.55])
plot(find(tp_cnt1_p<0.025),n_p_atx(find(tp_cnt1_p<0.025),2),'k.','markersize',30)
plot(find(tp_cnt1_n<0.025),n_p_atx(find(tp_cnt1_n<0.025),2),'k.','markersize',30)

plot(prctile(perm_n_n_atx(:,:,2),95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(perm_n_p_atx(:,:,2),95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,4); hold on
plot(n_p_dpz(:,2),'r-','linewidth',3)
plot(n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([1 14 -0.05 0.55])
plot(find(tp_cnt2_p<0.025),n_p_dpz(find(tp_cnt2_p<0.025),2),'k.','markersize',30)
plot(find(tp_cnt2_n<0.025),n_n_dpz(find(tp_cnt2_n<0.025),2),'k.','markersize',30)

plot(prctile(perm_n_n_dpz(:,:,2),95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(perm_n_p_dpz(:,:,2),95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,5); hold on
plot(d1_n,'b-','linewidth',3)
plot(d1_p,'r-','linewidth',3)
axis([0 14 -0.3 0.55])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(p_d1_p<0.025),d1_p(find(p_d1_p<0.025)),'k.','markersize',30)
plot(find(p_d1_n<0.025),d1_n(find(p_d1_n<0.025)),'k.','markersize',30)

subplot(3,2,6); hold on
plot(d2_n,'b-','linewidth',3)
plot(d2_p,'r-','linewidth',3)
axis([0 14 -0.3 0.55])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]');
plot(find(p_d2_p<0.025),d2_p(find(p_d2_p<0.025)),'k.','markersize',30)
plot(find(p_d2_n<0.025),d2_n(find(p_d2_n<0.025)),'k.','markersize',30)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp.pdf'));



%% MULTIPLE COMPARISONS CORRECTION (HAWELLEK ET AL., xxxx)

for ifreq = 1 : 13
    [~,~,idx_D_cnt1_n(:,ifreq)] = unique(perm_n_n_atx(:,ifreq,2));
    [~,~,idx_D_cnt1_p(:,ifreq)] = unique(perm_n_p_atx(:,ifreq,2));
    [~,~,idx_D_cnt2_n(:,ifreq)] = unique(perm_n_n_dpz(:,ifreq,2));
    [~,~,idx_D_cnt2_p(:,ifreq)] = unique(perm_n_p_dpz(:,ifreq,2));
    [~,~,idx_D_res1_n(:,ifreq)] = unique(perm_n_n_atx(:,ifreq,1));
    [~,~,idx_D_res1_p(:,ifreq)] = unique(perm_n_p_atx(:,ifreq,1));
    [~,~,idx_D_res2_n(:,ifreq)] = unique(perm_n_n_dpz(:,ifreq,1));
    [~,~,idx_D_res2_p(:,ifreq)] = unique(perm_n_p_dpz(:,ifreq,1));
    
    [~,~,idx_D_res1_all(:,ifreq)] = unique(perm_n_all_atx(:,ifreq,1));
    [~,~,idx_D_cnt1_all(:,ifreq)] = unique(perm_n_all_atx(:,ifreq,2));
    
    [~,~,idx_D_res2_all(:,ifreq)] = unique(perm_n_all_dpz(:,ifreq,1));
    [~,~,idx_D_cnt2_all(:,ifreq)] = unique(perm_n_all_dpz(:,ifreq,2));
    
   
end

idx_R_cnt1_p = max(idx_D_cnt1_p,[],2);
idx_R_cnt1_n = max(idx_D_cnt1_n,[],2);
idx_R_cnt2_p = max(idx_D_cnt2_p,[],2);
idx_R_cnt2_n = max(idx_D_cnt2_n,[],2);
idx_R_res1_n = max(idx_D_res1_n,[],2);
idx_R_res1_p = max(idx_D_res1_p,[],2);
idx_R_res2_n = max(idx_D_res2_n,[],2);
idx_R_res2_p = max(idx_D_res2_p,[],2);

idx_R_res1_all = max(idx_D_res1_all,[],2);
idx_R_res2_all = max(idx_D_res2_all,[],2);
idx_R_cnt1_all = max(idx_D_cnt1_all,[],2);
idx_R_cnt2_all = max(idx_D_res2_all,[],2);

for ifreq = 1 : 13
    for irank = 1 : length(idx_R_cnt1_p)
        irank
        
        tmp = perm_n_p_atx(find(idx_D_res1_p(:,ifreq) == idx_R_res1_p(irank)),ifreq,1);
        tperm_res1_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_res1_p(:,ifreq) == irank)),ifreq,1);
        tmp = perm_n_n_atx(find(idx_D_res1_n(:,ifreq) == idx_R_res1_n(irank)),ifreq,1);
        tperm_res1_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_res1_n(:,ifreq) == irank)),ifreq,1);
        
        tmp = perm_n_p_dpz(find(idx_D_res2_p(:,ifreq) == idx_R_res2_p(irank)),ifreq,1);
        tperm_res2_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_res2_p(:,ifreq) == irank)),ifreq,1);
        tmp = perm_n_n_dpz(find(idx_D_res2_n(:,ifreq) == idx_R_res2_n(irank)),ifreq,1);
        tperm_res2_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_res2_n(:,ifreq) == irank)),ifreq,1);
        
        tmp = perm_n_p_atx(find(idx_D_cnt1_p(:,ifreq) == idx_R_cnt1_p(irank)),ifreq,2);
        tperm_cnt1_p_corr(irank,ifreq)=perm_n_p_atx(unique(find(idx_D_cnt1_p(:,ifreq) == irank)),ifreq,2);
        tmp = perm_n_n_atx(find(idx_D_cnt1_n(:,ifreq) == idx_R_cnt1_n(irank)),ifreq,2);
        tperm_cnt1_n_corr(irank,ifreq)=perm_n_n_atx(unique(find(idx_D_cnt1_n(:,ifreq) == irank)),ifreq,2);
        
        tmp = perm_n_p_dpz(find(idx_D_cnt2_p(:,ifreq) == idx_R_cnt2_p(irank)),ifreq,2);
        tperm_cnt2_p_corr(irank,ifreq)=perm_n_p_dpz(unique(find(idx_D_cnt2_p(:,ifreq) == irank)),ifreq,2);
        tmp = perm_n_n_dpz(find(idx_D_cnt2_n(:,ifreq) == idx_R_cnt2_n(irank)),ifreq,2);
        tperm_cnt2_n_corr(irank,ifreq)=perm_n_n_dpz(unique(find(idx_D_cnt2_n(:,ifreq) == irank)),ifreq,2);
        
        
        % COMPUTE FOR ALL CONNECTIONS -------------------------------------
        % ATX, REST
        tmp = perm_n_all_atx(find(idx_D_res1_all(:,ifreq) == idx_R_res1_all(irank)),ifreq,1);
        tperm_res_atx_all_corr(irank,ifreq)=perm_n_all_atx(unique(find(idx_D_res1_all(:,ifreq) == irank)),ifreq,1);
        % DPZ, REST
        tmp = perm_n_all_dpz(find(idx_D_res2_all(:,ifreq) == idx_R_res2_all(irank)),ifreq,1);
        tperm_res_dpz_all_corr(irank,ifreq)=perm_n_all_atx(unique(find(idx_D_res2_all(:,ifreq) == irank)),ifreq,1);
%       % ATX, CNT
        tmp = perm_n_all_atx(find(idx_D_cnt1_all(:,ifreq) == idx_R_cnt1_all(irank)),ifreq,2);
        tperm_cnt_atx_all_corr(irank,ifreq)=perm_n_all_atx(unique(find(idx_D_cnt1_all(:,ifreq) == irank)),ifreq,2);
        % DPZ, CNT
        tmp = perm_n_all_dpz(find(idx_D_cnt2_all(:,ifreq) == idx_R_cnt2_all(irank)),ifreq,2);
        tperm_cnt_dpz_all_corr(irank,ifreq)=perm_n_all_atx(unique(find(idx_D_cnt2_all(:,ifreq) == irank)),ifreq,2);
%        
        % TASK VS REST ------------------------------------------------------
        [t_tvsr,~,~,s] = ttest(taskvsrest_perm(:,:,:,2,ifoi),taskvsrest_perm(:,:,:,1,ifoi),'dim',3,'alpha',alp);
        t_tvsr = t_tvsr.*sign(s.tstat); clear s   
      
        par.tperm_taskvsrest_p(kperm,ifoi) = nansum(nansum(t_tvsr>0))./8010;
        par.tperm_taskvsrest_n(kperm,ifoi) = nansum(nansum(t_tvsr<0))./8010;
      
      
        
    end
end

for ifreq = 1 : 13
    
    p_res1_p(ifreq) = 1-sum(n_p_atx(ifreq,1)>tperm_res1_p_corr(:,ifreq))/nperm;
    p_res1_n(ifreq) = 1-sum(n_n_atx(ifreq,1)>tperm_res1_n_corr(:,ifreq))/nperm;
    p_cnt1_p(ifreq) = 1-sum(n_p_atx(ifreq,2)>tperm_cnt1_p_corr(:,ifreq))/nperm;
    p_cnt1_n(ifreq) = 1-sum(n_n_atx(ifreq,2)>tperm_cnt1_n_corr(:,ifreq))/nperm;
    p_res2_p(ifreq) = 1-sum(n_p_dpz(ifreq,1)>tperm_res2_p_corr(:,ifreq))/nperm;
    p_res2_n(ifreq) = 1-sum(n_n_dpz(ifreq,1)>tperm_res2_n_corr(:,ifreq))/nperm;
    p_cnt2_p(ifreq) = 1-sum(n_p_dpz(ifreq,2)>tperm_cnt2_p_corr(:,ifreq))/nperm;
    p_cnt2_n(ifreq) = 1-sum(n_n_dpz(ifreq,2)>tperm_cnt2_n_corr(:,ifreq))/nperm;
    
    p_cnt_atx_all(ifreq) = 1-sum(atx(ifreq,2)>tperm_cnt_atx_all_corr(:,ifreq))/nperm;
    p_res_atx_all(ifreq) = 1-sum(atx(ifreq,1)>tperm_res_atx_all_corr(:,ifreq))/nperm;

    p_cnt_dpz_all(ifreq) = 1-sum(dpz(ifreq,2)>tperm_cnt_dpz_all_corr(:,ifreq))/nperm;
    p_res_dpz_all(ifreq) = 1-sum(dpz(ifreq,1)>tperm_res_dpz_all_corr(:,ifreq))/nperm;
    
end
%% ALTERNATIVE CORRECTION
% do not perform rank conversion, but only the nichols & holmes style

for ifreq = 1 : 13
  
  idx_R_cnt1_p(:,ifreq) = max(abs(perm_n_p_atx(:,:,2)),[],2);
  idx_R_cnt1_n(:,ifreq) = max(abs(perm_n_n_atx(:,:,2)),[],2);
  idx_R_cnt2_p(:,ifreq) = max(abs(perm_n_p_dpz(:,:,2)),[],2);
  idx_R_cnt2_n(:,ifreq) = max(abs(perm_n_n_dpz(:,:,2)),[],2);
  
  idx_R_res1_p(:,ifreq) = max(abs(perm_n_p_atx(:,:,1)),[],2);
  idx_R_res1_n(:,ifreq) = max(abs(perm_n_n_atx(:,:,1)),[],2);
  idx_R_res2_p(:,ifreq) = max(abs(perm_n_p_dpz(:,:,1)),[],2);
  idx_R_res2_n(:,ifreq) = max(abs(perm_n_n_dpz(:,:,1)),[],2);
  
  idx_R_context1_p(:,ifreq) = max(abs(perm_n_p_atx(:,:,1)-perm_n_p_atx(:,:,2)),[],2);
  idx_R_context1_n(:,ifreq) = max(abs(perm_n_n_atx(:,:,1)-perm_n_n_atx(:,:,2)),[],2);

  idx_R_context2_p(:,ifreq) = max(abs(perm_n_p_dpz(:,:,1)-perm_n_p_dpz(:,:,2)),[],2);
  idx_R_context2_n(:,ifreq) = max(abs(perm_n_n_dpz(:,:,1)-perm_n_n_dpz(:,:,2)),[],2);

  idx_R_res1_all(:,ifreq) = max(abs(perm_n_all_atx(:,:,1)),[],2);
  idx_R_res2_all(:,ifreq) = max(abs(perm_n_all_dpz(:,:,1)),[],2);
  idx_R_cnt1_all(:,ifreq) = max(abs(perm_n_all_atx(:,:,2)),[],2);
  idx_R_cnt2_all(:,ifreq) = max(abs(perm_n_all_dpz(:,:,2)),[],2);
  
  idx_R_doubledissociation(:,ifreq) = max(perm_doubledissociation,[],2);
  idx_R_taskvsrest_p(:,ifreq) = max(abs(taskvsrest_p_perm),[],2);
  idx_R_taskvsrest_n(:,ifreq) = max(abs(taskvsrest_n_perm),[],2);
  
  p_cnt1_p(ifreq) = 1-sum(abs(n_p_atx(ifreq,2))>abs(idx_R_cnt1_p(:,ifreq)))/nperm;
  p_cnt1_n(ifreq) = 1-sum(abs(n_n_atx(ifreq,2))>abs(idx_R_cnt1_n(:,ifreq)))/nperm;
  p_cnt2_p(ifreq) = 1-sum(abs(n_p_dpz(ifreq,2))>abs(idx_R_cnt2_p(:,ifreq)))/nperm;
  p_cnt2_n(ifreq) = 1-sum(abs(n_n_dpz(ifreq,2))>abs(idx_R_cnt2_n(:,ifreq)))/nperm;
  
  p_res1_p(ifreq) = 1-sum(abs(n_p_atx(ifreq,1))>abs(idx_R_res1_p(:,ifreq)))/nperm;
  p_res1_n(ifreq) = 1-sum(abs(n_n_atx(ifreq,1))>abs(idx_R_res1_n(:,ifreq)))/nperm;
  p_res2_p(ifreq) = 1-sum(abs(n_p_dpz(ifreq,1))>abs(idx_R_res2_p(:,ifreq)))/nperm;
  p_res2_n(ifreq) = 1-sum(abs(n_n_dpz(ifreq,1))>abs(idx_R_res2_n(:,ifreq)))/nperm;
    
  p_cnt_atx_all(ifreq) = 1-sum(abs(atx(ifreq,2))>abs(idx_R_cnt1_all(:,ifreq)))/nperm;
  p_res_atx_all(ifreq) = 1-sum(abs(atx(ifreq,1))>abs(idx_R_res1_all(:,ifreq)))/nperm;

  p_cnt_dpz_all(ifreq) = 1-sum(abs(dpz(ifreq,2))>idx_R_cnt2_all(:,ifreq))/nperm;
  p_res_dpz_all(ifreq) = 1-sum(abs(dpz(ifreq,1))>idx_R_res2_all(:,ifreq))/nperm;
    
  p_atx_context_p_corr(ifreq) = 1-sum([abs(n_p_atx(ifreq,1)-n_p_atx(ifreq,2))]>idx_R_context1_p(:,ifreq))/nperm;
  p_atx_context_n_corr(ifreq) = 1-sum([abs(n_n_atx(ifreq,1)-n_n_atx(ifreq,2))]>idx_R_context1_n(:,ifreq))/nperm;
    
  p_dpz_context_p_corr(ifreq) = 1-sum([abs(n_p_dpz(ifreq,1)-n_p_dpz(ifreq,2))]>idx_R_context2_p(:,ifreq))/nperm;
  p_dpz_context_n_corr(ifreq) = 1-sum([abs(n_n_dpz(ifreq,1)-n_n_dpz(ifreq,2))]>idx_R_context2_n(:,ifreq))/nperm;
    
  p_doubledissociation_corr(ifreq) = 1-sum(doubledissociation_emp(ifreq)>idx_R_doubledissociation(:,ifreq))/nperm;
  pval_taskvsrest_p_corr(ifreq) = 1-sum(abs(taskvsrest_p(ifreq))>idx_R_taskvsrest_p(:,ifreq))/nperm;
  pval_taskvsrest_n_corr(ifreq) = 1-sum(abs(taskvsrest_n(ifreq))>idx_R_taskvsrest_n(:,ifreq))/nperm;
  
  
end

%% PLOT CHANGES

ALPHA = 0.05;

figure;

subplot(3,2,1); hold on
plot(n_p_atx(:,1),'r-','linewidth',3)
plot(n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')
title('Rest')
axis([1 14 -0.05 0.55])
plot(find(p_res1_p<ALPHA),n_p_atx(find(p_res1_p<ALPHA),1),'k.','markersize',20)
plot(find(p_res1_n<ALPHA),n_n_atx(find(p_res1_n<ALPHA),1),'k.','markersize',20)

plot(prctile(idx_R_res1_n,95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(idx_R_res1_p,95),'linewidth',1,'color','r','linestyle',':')

subplot(3,2,2); hold on
plot(n_p_dpz(:,1),'r-','linewidth',3)
plot(n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([1 14 -0.05 0.55])
ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')
title('Rest')

plot(find(p_res2_p<ALPHA),n_p_dpz(find(p_res2_p<ALPHA),1),'k.','markersize',20)
plot(find(p_res2_n<ALPHA),n_n_dpz(find(p_res2_n<ALPHA),1),'k.','markersize',20)

plot(prctile(idx_R_res2_n,95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(idx_R_res2_p,95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,3); hold on
plot(n_p_atx(:,2),'r-','linewidth',3)
plot(n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')
title('Task')
axis([1 14 -0.05 0.55])
plot(find(p_cnt1_p<ALPHA),n_p_atx(find(p_cnt1_p<ALPHA),2),'k.','markersize',20)
plot(find(p_cnt1_n<ALPHA),n_p_atx(find(p_cnt1_n<ALPHA),2),'k.','markersize',20)

plot(prctile(idx_R_cnt1_n,95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(idx_R_cnt1_p,95),'linewidth',1,'color','r','linestyle',':')


subplot(3,2,4); hold on
plot(n_p_dpz(:,2),'r-','linewidth',3)
plot(n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4 0.5],'yticklabel',[0 10 20 30 40 50])
axis([1 14 -0.05 0.55])
title('Task')

ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')

plot(find(p_cnt2_p<ALPHA),n_p_dpz(find(p_cnt2_p<ALPHA),2),'k.','markersize',20)
plot(find(p_cnt2_n<ALPHA),n_n_dpz(find(p_cnt2_n<ALPHA),2),'k.','markersize',20)

plot(prctile(idx_R_cnt2_n,95),'linewidth',1,'color','b','linestyle',':')
plot(prctile(idx_R_cnt2_p,95),'linewidth',1,'color','r','linestyle',':')

% CONTEXT-DEPENDENCE ------------------------------------------------------
% ATOMOXETINE
subplot(3,2,5); hold on

dp = n_p_atx(:,1)-n_p_atx(:,2);
dn = n_n_atx(:,1)-n_n_atx(:,2);

plot(n_p_atx(:,1)-n_p_atx(:,2),'r-','linewidth',3)
plot(n_n_atx(:,1)-n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-0.5 -0.25 0 0.25 0.5],'yticklabel',[-50 -25 0 25 50])
ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')
title('Context')
axis([1 14 -0.6 0.6])
plot(find(p_atx_context_p_corr<ALPHA),dp(find(p_atx_context_p_corr<ALPHA)),'k.','markersize',20)
plot(find(p_atx_context_n_corr<ALPHA),dn(find(p_atx_context_n_corr<ALPHA)),'k.','markersize',20)

% plot(prctile(idx_R_cnt1_n,95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(idx_R_cnt1_p,95),'linewidth',1,'color','r','linestyle',':')

% DONEPEZIL
subplot(3,2,6); hold on
dp = n_p_dpz(:,1)-n_p_dpz(:,2);
dn = n_n_dpz(:,1)-n_n_dpz(:,2);

plot(n_p_dpz(:,1)-n_p_dpz(:,2),'r-','linewidth',3)
plot(n_n_dpz(:,1)-n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-0.5 -0.25 0 0.25 0.5],'yticklabel',[-50 -25 0 25 50])
ylabel('Altered corr. [%]'); xlabel('Frequency [Hz]')
title('Context')
axis([1 14 -0.6 0.6])
plot(find(p_dpz_context_p_corr<ALPHA),dp(find(p_dpz_context_p_corr<ALPHA)),'k.','markersize',20)
plot(find(p_dpz_context_n_corr<ALPHA),dn(find(p_dpz_context_n_corr<ALPHA)),'k.','markersize',20)

% plot(prctile(idx_R_cnt1_n,95),'linewidth',1,'color','b','linestyle',':')
% plot(prctile(idx_R_cnt1_p,95),'linewidth',1,'color','r','linestyle',':')


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_alteredcorrelations_bothdir_corrected.pdf'))


figure;
% DONEPEZIL
subplot(3,2,6); hold on

plot(doubledissociation_emp,'color','[0.7 0.7 0.7]','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-0.5 -0.25 0 0.25 0.5],'yticklabel',[-50 -25 0 25 50])
ylabel('Double difference [%]'); xlabel('Frequency [Hz]')
title('Double dissociation')
axis([1 14 -0.75 0.75])
plot(find(p_doubledissociation_corr<ALPHA),doubledissociation_emp(find(p_doubledissociation_corr<ALPHA)),'k.','markersize',20)

% plot(prctile(idx_R_cnt1_n,95),'linewidth',1,'color','b','linestyle',':')
% plo
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_doubledissociation.pdf'))

%% PLOT TASK VS REST, FC MATRIX AND SPECTRUM

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip'); cmap = cmap(end:-1:1,:); 

figure;

subplot(3,2,3); hold on
d = nanmean(cleandat(:,:,:,1,2,7),3)./nanmean(cleandat(:,:,:,1,1,7),3);
d=triu(d,1);d(d==0)=1;

d = rot90(d);

imagesc(d,[0.70 1.3]);
colormap(cmap); axis square; box off; axis off

subplot(3,2,5); hold on

d = nanmean(cleandat(:,:,:,1,2,6),3)./nanmean(cleandat(:,:,:,1,1,6),3);
d=triu(d,1);d(d==0)=1;

d = rot90(d);

imagesc(d,[0.70 1.3]);
colormap(cmap); axis square; box off; axis off


subplot(3,2,6); hold on

plot(taskvsrest_p,'color','[0.7 0.7 0.7]','linewidth',3)
plot(taskvsrest_n,'color','[0.7 0.7 0.7]','linewidth',3)

set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[-0.5 -0.25 0 0.25 0.5],'yticklabel',[-50 -25 0 25 50])
ylabel(sprintf('Fraction of \n altered correlations [%%]')); xlabel('Frequency [Hz]')
title('Double dissociation')
axis([1 14 -0.05 0.5])
plot(find(pval_taskvsrest_n_corr<ALPHA),taskvsrest_n(find(pval_taskvsrest_n_corr<ALPHA)),'k.','markersize',20)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_taskvsrest.pdf'))