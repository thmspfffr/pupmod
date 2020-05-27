% %% pupmod_all_dimensionality
% COMPUTES DIMENSIONALITY OF FC MATRICES.

% Last changed: 13-11-2018

clear

v = 3;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

fprintf('Loading grid...\n')
grid  = select_chans(sa_meg_template.grid_cortex3000,400); fprintf('Loading grid... Done\n')
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

fc = pupmod_loadcov(v,SUBJLIST,0);

if v == 20
  tmp = tp_create_grid('vtpm');
  clim = [-0.025 0.025];
  idx = 1:46;
  iismember(idx,[21 22 23 44 45 46]);
  reg = tmp.tissuelabel_4mm(idx);
end


%% COMPUTE DIMENSIONALITY
dim = zeros(28,3,2,21);
for isubj = 1:28
  isubj
  cnt = 0;
  for m = 1 : 3
    for ifoi=1:21
      for cond = 2 : 2
%         cnt=cnt+1
        % compute according to ito et al. (2019) biorxiv
        [a,b,eig]=pca(nanmean(fc(:,:,isubj,m,cond,ifoi,:),7));
        dim(isubj,m,cond,ifoi)=(sum(eig)^2)/(sum(eig.^2));
        
      end
    end
  end
end

%% PLOT

figure; set(gcf,'color','w');
ipharm = 1;
subplot(3,2,1); hold on; title('Placebo')
plot(squeeze(nanmean(nanmean(dim(:,ipharm,1:2,:)),2))')
legend({'Rest';'Task'}); legend boxoff
[h1,p1]=ttest(squeeze(nanmean(dim(:,ipharm,2,:),2)),squeeze(nanmean(dim(:,ipharm,1,:),2)),'dim',1);
k=plot(find(h1),ones(sum(h1)),'.','markersize',5,'color','k'),
% set(get(get(k,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

ipharm = 2;
subplot(3,2,3); hold on; title('Atomoxetine')
plot(squeeze(nanmean(nanmean(dim(:,ipharm,1:2,:)),2))')
[h2,p2]=ttest(squeeze(nanmean(dim(:,ipharm,2,:),2)),squeeze(nanmean(dim(:,ipharm,1,:),2)),'dim',1);
plot(find(h2),ones(sum(h2)),'.','markersize',5,'color','k'),
tp_editplots; xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

ipharm = 3;
subplot(3,2,5); hold on; title('Donepezil')
plot(squeeze(nanmean(nanmean(dim(:,ipharm,1:2,:)),2))')
[h3,p3]=ttest(squeeze(nanmean(dim(:,ipharm,2,:),2)),squeeze(nanmean(dim(:,ipharm,1,:),2)),'dim',1);
plot(find(h3),ones(sum(h3)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

subplot(3,2,2); hold on; title('Average')
plot(squeeze(nanmean(nanmean(dim(:,1:3,1:2,:)),2))')
[h_avg,p_avg]=ttest(squeeze(nanmean(dim(:,1:3,2,:),2)),squeeze(nanmean(dim(:,1:3,1,:),2)),'dim',1);
plot(find(h_avg),ones(sum(h_avg)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

subplot(3,2,4); hold on; title('Atx vs Pbo (Task)')
plot(squeeze(nanmean(dim(:,[1 2],2,:)))')
[h_atx,p_atx]=ttest(squeeze(nanmean(dim(:,2,2,:),2)),squeeze(nanmean(dim(:,1,2,:),2)),'dim',1);
plot(find(h_atx),ones(sum(h_atx)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

subplot(3,2,6); hold on; title('Dpz vs Pbo (Rest)')
plot(squeeze(nanmean(dim(:,[1 3],1,:)))')
[h_dpz,p_dpz]=ttest(squeeze(nanmean(dim(:,3,1,:),2)),squeeze(nanmean(dim(:,1,1,:),2)),'dim',1);
plot(find(h_dpz),ones(sum(h_dpz)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 0 15])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_dimensionality_v%d.pdf',v))

%%
behav_pooled = nanmean(behav_pooled,3);
behav = nanmean(behav_cnt,3);
d_behav = nanmean(behav_cnt(2,:,:),3)-nanmean(behav_cnt(1,:,:),3);

% d_behav = behav_pooled(2,:)-behav_pooled(1,:,:);

d_dim2 = squeeze(dim(~isnan(d_behav),2,2,:)-dim(~isnan(d_behav),1,2,:));
d_dim1 = squeeze(dim(~isnan(d_behav),2,1,:)-dim(~isnan(d_behav),1,1,:));

[rr1,p1]=corr(d_dim1,d_behav(~isnan(d_behav))')
[rr2,p2]=corr(d_dim2,d_behav(~isnan(d_behav))')

perm = 1;


nperm= 20000;
   all_idx1 = randi(2,[28,nperm]);

for iperm =1  :nperm
iperm
    idx1 = all_idx1(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(i,1) = behav(idx1(i),i);
      permdat(i,2) = behav(idx2(i),i);
    end
    
  d_behav = permdat(:,2)-permdat(:,1);

  [r1(:,iperm)]=corr(d_dim1,d_behav);  
  [r2(:,iperm)]=corr(d_dim2,d_behav);

  
  
end

%% PLOT MEAN FC AS COMPARISON
fc1 = squeeze(nanmean(nanmean(nanmean(fc,1),2),7));

figure; set(gcf,'color','w');
ipharm = 1;
subplot(3,2,1); hold on; title('Placebo')
plot(squeeze(nanmean(nanmean(fc1(:,ipharm,1:2,:)),2))')
legend({'Rest';'Task'}); legend boxoff
[h1,p1]=ttest(squeeze(nanmean(fc1(:,ipharm,2,:),2)),squeeze(nanmean(fc1(:,ipharm,1,:),2)),'dim',1);
k=plot(find(h1),0*ones(sum(h1)),'.','markersize',5,'color','k'),
set(get(get(k,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

ipharm = 2;
subplot(3,2,3); hold on; title('Atomoxetine')
plot(squeeze(nanmean(nanmean(fc1(:,ipharm,1:2,:)),2))')
[h2,p2]=ttest(squeeze(nanmean(fc1(:,ipharm,2,:),2)),squeeze(nanmean(fc1(:,ipharm,1,:),2)),'dim',1);
plot(find(h2),0*ones(sum(h2)),'.','markersize',5,'color','k'),
tp_editplots; xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

ipharm = 3;
subplot(3,2,5); hold on; title('Donepezil')
plot(squeeze(nanmean(nanmean(fc1(:,ipharm,1:2,:)),2))')
[h3,p3]=ttest(squeeze(nanmean(fc1(:,ipharm,2,:),2)),squeeze(nanmean(fc1(:,ipharm,1,:),2)),'dim',1);
plot(find(h3),0*ones(sum(h3)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

subplot(3,2,2); hold on; title('Average')
plot(squeeze(nanmean(nanmean(fc1(:,1:3,1:2,:)),2))')
[h_avg,p_avg]=ttest(squeeze(nanmean(fc1(:,1:3,2,:),2)),squeeze(nanmean(fc1(:,1:3,1,:),2)),'dim',1);
plot(find(h_avg),0*ones(sum(h_avg)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

subplot(3,2,4); hold on; title('Atx vs Pbo (Task)')
plot(squeeze(nanmean(fc1(:,[1 2],2,:)))')
[h_atx,p_atx]=ttest(squeeze(nanmean(fc1(:,2,2,:),2)),squeeze(nanmean(fc1(:,1,2,:),2)),'dim',1);
plot(find(h_atx),0*ones(sum(h_atx)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

subplot(3,2,6); hold on; title('Dpz vs Pbo (Rest)')
plot(squeeze(nanmean(fc1(:,[1 3],1,:)))')
[h_dpz,p_dpz]=ttest(squeeze(nanmean(fc1(:,3,1,:),2)),squeeze(nanmean(fc1(:,1,1,:),2)),'dim',1);
plot(find(h_dpz),0*ones(sum(h_dpz)),'.','markersize',5,'color','k'),
tp_editplots;  xlabel('Frequency [Hz]'); ylabel('W')
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21 25],'xticklabel',[2 4 8 16 32 64 128])
axis([0 25 -0.01 .15])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_dimensionality_fc_v%d.pdf',v))

