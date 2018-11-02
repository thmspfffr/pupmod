%% PLOT NUMBER OF ALTERED CORRELATIONS, INCLUDING STATISTICS
% --------------------------
% This script obtains the empircal number of altered correlations, by
% calling pupmod_compute_altered_correlations.m and obtains a corrected 
% p-values from a permutation distribution (computed in
% pupmod_src_powcorr_permtest.m). The actual p-values are obtained calling
% the function pupmod_all_powcorr_getstatistics.m)
% --------------------------

clear

% -------------
% version of cleaned data: 
% v1: AAL, v2: 400 vertices (cortex)
% -------------
v = 12;
% -------------

% load  data
cd ~/pupmod/matlab/
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
%%

para = [];
para.nfreq = 13;
para.alpha = 0.05;

emp = pupmod_compute_altered_correlations(cleandat,para);

%% PLOT: No stats
% ------------------

linewidth = 2;

figure; set(gcf,'color','w');

subplot(3,2,1); hold on

plot(emp.n_p_atx(:,1),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
ylabel('Frac. of altered corr. [%]')

subplot(3,2,3); hold on

plot(emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
ylabel('Frac. of altered corr. [%]')

subplot(3,2,2); hold on

plot(emp.n_p_dpz(:,1),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,4); hold on

plot(emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -0.05 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

subplot(3,2,5); hold on

plot(emp.n_p_atx(:,1)-emp.n_p_atx(:,2),'r','linewidth',linewidth)
plot(emp.n_n_atx(:,1)-emp.n_n_atx(:,2),'b','linewidth',linewidth)

axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')

subplot(3,2,6); hold on

plot(emp.n_p_dpz(:,1)-emp.n_p_dpz(:,2),'r','linewidth',linewidth)
plot(emp.n_n_dpz(:,1)-emp.n_p_dpz(:,2),'b','linewidth',linewidth)

axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel('Frac. of altered corr. [%]')

%%
if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
  % Settings:
  para = [];
  para.nfreq = [1:13]; % freqs 1 - 13
  para.alpha = 0.05; % alpha for adjacency
  para.ver = v;
  para.nperm = 50000;
  para.nsubs = 100;
  para.type = 'global';
  para.cond = 'atx';
  [outp_atx, emp_atx] = pupmod_src_powcorr_getstatistics(para);
  para.cond = 'dpz';
  [outp_dpz, emp_dpz] = pupmod_src_powcorr_getstatistics(para);
  
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v),'outp_atx','outp_dpz','emp_atx','emp_dpz')
else
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_alteredcorr_v%d.mat',v))
end
%% PLOT P-VALUES
figure;  set(gcf,'color','white')

subplot(3,2,1); hold on
plot(-log10(outp.p_res1_n),'b-','linewidth',3)
plot(-log10(outp.p_res1_p),'r-','linewidth',3)

line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Rest')

subplot(3,2,2); hold on
plot(-log10(outp.p_res2_n),'b-','linewidth',3)
plot(-log10(outp.p_res2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

subplot(3,2,3); hold on
plot(-log10(outp.p_cnt1_n),'b-','linewidth',3)
plot(-log10(outp.p_cnt1_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
title('Task')

subplot(3,2,4); hold on
plot(-log10(outp.p_cnt2_n),'b-','linewidth',3)
plot(-log10(outp.p_cnt2_p),'r-','linewidth',3)
line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
axis([0 14 0 4])
set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')

% subplot(3,2,5); hold on
% plot(-log10(p_d1_n),'b-','linewidth',3)
% plot(-log10(p_d1_p),'r-','linewidth',3)
% line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
% axis([0 14 0 4])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
% 
% subplot(3,2,6); hold on
% plot(-log10(p_d2_n),'b-','linewidth',3)
% plot(-log10(p_d2_p),'r-','linewidth',3)
% line([0 14],[-log10(0.025) -log10(0.025)],'linestyle','--','color','k')
% axis([0 14 0 4])
% set(gca,'tickdir','out','ytick',[0 1 2 3],'yticklabel',[0 0.1 0.01 0.001])
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% xlabel('Carrier frequency [Hz]'); ylabel('P-Value (uncorrected)')
% % tp_editplots
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_powcorr_taskrestcomp_pval.pdf'));

%%

alpha1 = 0.025;
alpha2 = 0.01;
alpha3 = 0.001;

figure; set(gcf,'color','w')

subplot(3,2,1); hold on
plot(emp_atx.n_p_atx(:,1),'r-','linewidth',3)
plot(emp_atx.n_n_atx(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.5],'yticklabel',[0 25 50])
ylabel('Altered corr. [%]')
title('Rest')
axis([0 14 -0.05 0.55])
plot(find(outp_atx.p_res1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_res1_p<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_res1_n<alpha3),1),'ro','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(3,2,2); hold on
plot(emp_atx.n_p_dpz(:,1),'r-','linewidth',3)
plot(emp_atx.n_n_dpz(:,1),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.5],'yticklabel',[0 25 50])
axis([0 14 -0.05 0.55])
plot(find(outp_atx.p_res2_p<alpha1),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_n<alpha1),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha1),1),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_res2_p<alpha2),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_n<alpha2),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha2),1),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_res2_p<alpha3),emp_atx.n_p_dpz(find(outp_atx.p_res2_p<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_res2_n<alpha3),emp_atx.n_n_dpz(find(outp_atx.p_res2_n<alpha3),1),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(3,2,3); hold on
plot(emp_atx.n_p_atx(:,2),'r-','linewidth',3)
plot(emp_atx.n_n_atx(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.5],'yticklabel',[0 25 50])
ylabel('Altered corr. [%]')
title('Task')
axis([0 14 -0.05 0.55])
plot(find(outp_atx.p_cnt1_p<alpha1),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_n<alpha1),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt1_p<alpha2),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_n<alpha2),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt1_p<alpha3),emp_atx.n_p_atx(find(outp_atx.p_cnt1_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt1_n<alpha3),emp_atx.n_n_atx(find(outp_atx.p_cnt1_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(3,2,4); hold on
plot(emp_atx.n_p_dpz(:,2),'r-','linewidth',3)
plot(emp_atx.n_n_dpz(:,2),'b-','linewidth',3)
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
set(gca,'tickdir','out','ytick',[0 0.25 0.5],'yticklabel',[0 25 50])
axis([0 14 -0.05 0.55])
plot(find(outp_atx.p_cnt2_p<alpha1),emp_atx.n_p_dpz(find(outp_atx.p_cnt2_p<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_n<alpha1),emp_atx.n_n_dpz(find(outp_atx.p_cnt2_n<alpha1),2),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_cnt2_p<alpha2),emp_atx.n_p_dpz(find(outp_atx.p_cnt2_p<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_n<alpha2),emp_atx.n_n_dpz(find(outp_atx.p_cnt2_n<alpha2),2),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_cnt2_p<alpha3),emp_atx.n_p_dpz(find(outp_atx.p_cnt2_p<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_cnt2_n<alpha3),emp_atx.n_n_dpz(find(outp_atx.p_cnt2_n<alpha3),2),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(3,2,5); hold on
plot(emp_atx.n_n_context_atx,'b-','linewidth',3)
plot(emp_atx.n_p_context_atx,'r-','linewidth',3)
axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','ytick',[-0.5 0 0.5],'yticklabel',[-50 0 50])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); ylabel('Difference')
plot(find(outp_atx.p_context1_p<alpha1),emp_atx.n_p_context_atx(find(outp.p_context1_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_n<alpha1),emp_atx.n_n_context_atx(find(outp.p_context1_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context1_p<alpha2),emp_atx.n_p_context_atx(find(outp.p_context1_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_n<alpha2),emp_atx.n_n_context_atx(find(outp.p_context1_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context1_p<alpha3),emp_atx.n_p_context_atx(find(outp.p_context1_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context1_n<alpha3),emp_atx.n_n_context_atx(find(outp.p_context1_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

subplot(3,2,6); hold on
plot(emp_atx.n_n_context_dpz,'b-','linewidth',3)
plot(emp_atx.n_p_context_dpz,'r-','linewidth',3)
axis([0 14 -0.6 0.6])
set(gca,'tickdir','out','ytick',[-0.5 0 0.5],'yticklabel',[-50 0 50])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Carrier frequency [Hz]'); 
plot(find(outp_atx.p_context2_p<alpha1),emp_atx.n_p_context_dpz(find(outp.p_context2_p<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_n<alpha1),emp_atx.n_n_context_dpz(find(outp.p_context2_n<alpha1)),'ko','markersize',7,'markerfacecolor','k')
plot(find(outp_atx.p_context2_p<alpha2),emp_atx.n_p_context_dpz(find(outp.p_context2_p<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_n<alpha2),emp_atx.n_n_context_dpz(find(outp.p_context2_n<alpha2)),'ko','markersize',7,'markerfacecolor','w')
plot(find(outp_atx.p_context2_p<alpha3),emp_atx.n_p_context_dpz(find(outp.p_context2_p<alpha3)),'ko','markersize',7,'markerfacecolor','m')
plot(find(outp_atx.p_context2_n<alpha3),emp_atx.n_n_context_dpz(find(outp.p_context2_n<alpha3)),'ko','markersize',7,'markerfacecolor','m')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_powcorr_drugeffects_lineplots.pdf'));

%% CREATE CIRCULAR PLOT
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',1));


mask = logical(tril(ones(76,76),-1));
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

nperm = 1000;

for iperm = 1 : nperm
    
    % within subjects permutation test 
    disp(sprintf('Perm #%d',iperm));
    
    idx1 = randi(2,[28 1]);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(:,:,i,1) = cleandat(:,:,i,idx1(i),2,6);
      permdat(:,:,i,2) = cleandat(:,:,i,idx2(i),2,6);     
    end
    
    [~,~,~,s]=ttest(permdat(include_bcn,include_bcn,:,2),permdat(include_bcn,include_bcn,:,1),'dim',3);
    maxt_atx(iperm) = max(abs(s.tstat(mask)));
        
    idx1 = randi(2,[28 1]); 
    idx2 = 3-idx1; 
    idx1(idx1==2)=3;
    idx2(idx2==2)=3;
    
    for i = 1 : length(idx1)
      permdat(:,:,i,1) = cleandat(:,:,i,idx1(i),1,7);
      permdat(:,:,i,2) = cleandat(:,:,i,idx2(i),1,7);     
    end
    
    [~,~,~,s]=ttest(permdat(include_bcn,include_bcn,:,2),permdat(include_bcn,include_bcn,:,1),'dim',3);
    maxt_dpz(iperm) = max(abs(s.tstat(mask)));
  
end
%%
[h,p,~,s] = ttest(cleandat(include_bcn,include_bcn,:,2,2,6),cleandat(include_bcn,include_bcn,:,1,2,6),'dim',3);
mean_change = nanmean(cleandat(include_bcn,include_bcn,:,2,2,6)-cleandat(include_bcn,include_bcn,:,1,2,6),3);
conn = s.tstat>prctile(maxt_atx,95);
conn = p < 0.0001;


labels = aal_labels_bcn; 
labels = labels(include_bcn);
labels(sum(conn)==0) = {''}


subplot(1,2,1);
circularGraph(mean_change.*conn,'label',labels)


[h,p,~,s] = ttest(cleandat(include_bcn,include_bcn,:,3,1,7),cleandat(include_bcn,include_bcn,:,1,1,7),'dim',3);
mean_change = nanmean(cleandat(include_bcn,include_bcn,:,3,1,7)-cleandat(include_bcn,include_bcn,:,1,1,7),3);
conn = abs(s.tstat)>prctile(maxt_dpz,95);

labels = aal_labels_bcn; 
labels = labels(include_bcn);
labels(sum(conn)==0) = {''}

subplot(1,2,2);
circularGraph(mean_change.*conn,'label',labels)
