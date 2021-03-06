%% ATX: Correlations of drug effects during Rest and Task
% compare drug effect during rest with drug effect during task

mask = logical(triu(ones(400,400),1));
fc = pupmod_loadpowcorr(23,SUBJLIST,0);
% fc = mod_loadpowcorr(25,SUBJLIST,1);

% ATX Drug effects
%%
for isubj = 1 : 28
  
  % Rest
  tmp = nanmean(nanmean(fc(:,:,isubj,:,:,11,:),7),6);
  tmp_rest = tmp(:,:,:,2,1)-tmp(:,:,:,1,1);
  tmp_rest = tmp_rest(mask);
  tmp_task = tmp(:,:,:,3,2)-tmp(:,:,:,1,2);
  tmp_task = tmp_task(mask);
  
  r_atx(isubj) = corr(tmp_rest,tmp_task);
  
  % Task
  tmp_rest = tmp(:,:,:,2,1)-tmp(:,:,:,1,1);
  tmp_rest = tmp_rest(mask);
  tmp_task = tmp(:,:,:,2,2)-tmp(:,:,:,1,2);
  tmp_task = tmp_task(mask);
  
  r_dpz(isubj) = corr(tmp_rest,tmp_task);
  
end
%%
% Correlation of mean effects (across subjects)
% --------------------------
tmp = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,[1 2],[1 2],10:12,:),2),1),7),6));

% Rest
tmp_rest_atx = squeeze(tmp(:,2,1)-tmp(:,1,1));
% Task
tmp_task_atx = squeeze(tmp(:,2,2)-tmp(:,1,2));

[r_atx_mean, p_atx_mean] = corr(tmp_rest_atx,tmp_task_atx);

tmp = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,[1 3],[1 2],13,:),2),1),7));
% Rest
tmp_rest_dpz  = squeeze(tmp(:,2,1)-tmp(:,1,1));
% Task
tmp_task_dpz = squeeze(tmp(:,2,2)-tmp(:,1,2));

[r_dpz_mean, p_dpz_mean] = corr(tmp_rest_dpz,tmp_task_dpz);


%%
figure; set(gcf,'color','w');

subplot(2,3,1); hold on

scatter(ones(28,1)-(rand(28,1)-0.5)/10,r_atx,'bo')
scatter(1.2,mean(r_atx),1000,'b.')

scatter(1.5*ones(28,1)-(rand(28,1)-0.5)/10,r_dpz,'ro')
scatter(1.7,mean(r_dpz),1000,'r.')

line([0.5 2], [0 0],'linestyle',':','color','k')

axis([0.5 2 -0.62 0.62]); tp_editplots;
ylabel('Correlation (FC mats)'); 

subplot(2,3,2); hold on
scatter(tmp_rest_atx,tmp_task_atx)
axis([-0.05 0.05 -0.05 0.05]); 
lsline; tp_editplots

subplot(2,3,3); hold on
scatter(tmp_rest_dpz,tmp_task_dpz)
axis([-0.05 0.05 -0.05 0.05]); 
lsline; tp_editplots

%% CORRELATION OF TASK AND DRUG EFFECTS
rand_fc = rand(400,400,28,2,2);

for isubj = 1 : 28
  
  % Task effect
  tmp = fc(:,:,isubj,[1],[1 2],6);
  tmp_task = tmp(:,:,:,1,2)-tmp(:,:,:,1,1);
  tmp_task = tmp_task(mask);
  
  % Atx
  tmp = fc(:,:,isubj,[1 2],[2 2],6);
  tmp_atx = tmp(:,:,:,2,2)-tmp(:,:,:,1,2);
  tmp_atx = tmp_atx(mask);
  
  r_atx(isubj) = corr(tmp_task,tmp_atx);
      
  % Task effect
  tmp = fc(:,:,isubj,[1],[1 2],7);
  tmp_task = tmp(:,:,:,1,2)-tmp(:,:,:,1,1);
  tmp_task = tmp_task(mask);
  
  % Atx
  tmp = fc(:,:,isubj,[1 3],[2 2],7);
  tmp_dpz = tmp(:,:,:,2,2)-tmp(:,:,:,1,2);
  tmp_dpz = tmp_dpz(mask);
  
  r_dpz(isubj) = corr(tmp_task,tmp_dpz);
  
end

%% CORRELATION ACROSS BLOCKS

% atx effect block 1
for ifoi = 1 : 25
  for isubj = 1 : 28

  % ATX effect
    tmp = fc(:,:,isubj,[1 2],2,ifoi,[1 2]);
    tmp1 = tmp(:,:,:,2,1,1)-tmp(:,:,:,1,1,1);
    tmp1 = tmp1(mask);

    tmp2 = tmp(:,:,:,2,1,2)-tmp(:,:,:,1,1,2);
    tmp2 = tmp2(mask);

    r12_atx(isubj,ifoi) = corr(tmp1,tmp2);

    tmp = fc(:,:,isubj,[1 3],1,ifoi,[1 2]);
    tmp1 = tmp(:,:,:,2,1,1)-tmp(:,:,:,1,1,1);
    tmp1 = tmp1(mask);

    tmp2 = tmp(:,:,:,2,1,2)-tmp(:,:,:,1,1,2);
    tmp2 = tmp2(mask);

    r12_dpz(isubj,ifoi) = corr(tmp1,tmp2);


  end
end

%% BASLINE EFFECTS

baseline = squeeze(nanmean(fc(:,:,isubj,[1],2,13,[1 2]),6));

for isubj = 1 : 28
  
  tmp = squeeze(nanmean(fc(:,:,isubj,[1 2],2,13,[1 2]),6));
  
  contr1 = tmp(:,:,2,1)-tmp(:,:,1,1);
  contr2 = tmp(:,:,2,2)-tmp(:,:,1,2);
  
  bl1 = baseline(:,:,1);
  bl2 = baseline(:,:,2);
  
  r1(isubj) = corr(contr1(mask),bl2(mask)); 
  r2(isubj) = corr(contr2(mask),bl1(mask));
  
  bl_mean = nanmean(baseline,3);
  tmp = nanmean(tmp,4);
  contr_mean = tmp(:,:,2)-tmp(:,:,1);
  r_mean(isubj) = corr(contr_mean(mask),bl_mean(mask));
  
end


%% ARE EFFECTS CORRELATED ACROSS SUBJECTS??
clear tmp
for isubj = 1 : 28
  tmp_fc = nanmean(fc(:,:,isubj,2,2,12,:),7)-nanmean(fc(:,:,isubj,1,2,12,:),7);
  tmp(:,isubj) = tmp_fc(mask);
end


figure; set(gcf,'color','w');
imagesc(corr(tmp),[-0.2 0.2]); axis square; set(gca,'ydir','normal')
xlabel('Subjects 1-28'); ylabel('Subjects 1-28'); tp_editplots

%% MEAN EFFECT CORRELATED ACROSS BLOCKS?
% ----------------------

figure; set(gcf,'color','w');

 ifoi = 10:12
% ATX
clear tmp 
tmp_fc1 = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,2,2,ifoi,1)-fc(:,:,:,1,2,ifoi,1),1),2),6));
tmp_fc2 = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,2,2,ifoi,2)-fc(:,:,:,1,2,ifoi,2),1),2),6));

[r_atx,p_atx]=corr(tmp_fc1(~isnan(tmp_fc1)),tmp_fc2(~isnan(tmp_fc1)))
% 
% subplot(1,2,1); hold on
% scatter(tmp_fc1(~isnan(tmp_fc1)),tmp_fc2(~isnan(tmp_fc1)));
% line([0 0],[-0.1 0.1],'color','k','linestyle',':')
% line([-0.1 0.1],[0 0],'color','k','linestyle',':')
% axis([-0.1 0.1 -0.1 0.1]); lsline; axis square; tp_editplots
% title(sprintf('ATX: r = %.3f | p = %.3f',r_atx,p_atx))

 ifoi = 13:14

% DPZ
clear tmp
tmp_fc1 = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,3,1,ifoi,1)-fc(:,:,:,1,1,ifoi,1),1),2),6));
tmp_fc2 = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,3,1,ifoi,2)-fc(:,:,:,1,1,ifoi,2),1),2),6));

[r_dpz,p_dpz]=corr(tmp_fc1(~isnan(tmp_fc1)),tmp_fc2(~isnan(tmp_fc1)));

% 
% subplot(1,2,2); hold on
% scatter(tmp_fc1(~isnan(tmp_fc1)),tmp_fc2(~isnan(tmp_fc1)),'ro')
% line([0 0],[-0.1 0.1],'color','k','linestyle',':')
% line([-0.1 0.1],[0 0],'color','k','linestyle',':')
% axis([-0.1 0.1 -0.1 0.1]); lsline; axis square; tp_editplots
% title(sprintf('DPZ: r = %.3f | p = %.3f',r_dpz,p_dpz))

% ------------------

%% MEAN FC EFFECT FOR THE TWO BLOCKS SEPARATELY
% --------------------
% ++++++++++++++++++++
% --------------------

atx_task = squeeze(nanmean(nanmean(fc(:,:,:,2,2,6,:))))-squeeze(nanmean(nanmean(fc(:,:,:,1,2,6,:))));
dpz_rest = squeeze(nanmean(nanmean(fc(:,:,:,3,1,7,:))))-squeeze(nanmean(nanmean(fc(:,:,:,1,1,7,:))));
atx_rest = squeeze(nanmean(nanmean(fc(:,:,:,2,1,6,:))))-squeeze(nanmean(nanmean(fc(:,:,:,1,1,6,:))));
dpz_task = squeeze(nanmean(nanmean(fc(:,:,:,3,2,7,:))))-squeeze(nanmean(nanmean(fc(:,:,:,1,2,7,:))));

figure; set(gcf,'color','w');
pos = ones(28,1)-((rand(28,1)-0.5)/10);

% --------------------
% ATOMOXETINE RESULTS DURING REST
% --------------------
subplot(2,2,1); hold on

for isubj =1 : 28
  c(isubj,:) = rand(1,3);
  scatter(pos(isubj),atx_rest(isubj,1),'w','markerfacecolor',c(isubj,:))
  scatter(1+pos(isubj),atx_rest(isubj,2),'w','markerfacecolor',c(isubj,:))
  line([pos(isubj) 1+pos(isubj)],[atx_rest(isubj,1) atx_rest(isubj,2)],'color',c(isubj,:))
end

line([-0.5 2.5],[0 0],'linestyle',':','color','k')
axis([0.5 2.5 -0.1 0.1]); tp_editplots

[h,p]=ttest(atx_rest);
text(1,0.075,sprintf('p = %.2f',p(1)),'horizontalalignment','center')
text(2,0.075,sprintf('p = %.2f',p(2)),'horizontalalignment','center')
title('Mean FC (Atx vs. Pbo): Rest'); 
set(gca,'xtick',[1 2],'xticklabel',{'Block1';'Block2'})
% --------------------
% DONEPEZIL RESULTS DURING REST
% --------------------
subplot(2,2,2); hold on

for isubj = 1 : 28
  scatter(pos(isubj),dpz_rest(isubj,1),'w','markerfacecolor',c(isubj,:))
  scatter(1+pos(isubj),dpz_rest(isubj,2),'w','markerfacecolor',c(isubj,:))
  line([pos(isubj) 1+pos(isubj)],[dpz_rest(isubj,1) dpz_rest(isubj,2)],'color',c(isubj,:))
end

line([-0.5 2.5],[0 0],'linestyle',':','color','k')
axis([0.5 2.5 -0.1 0.1]); tp_editplots

[h,p]=ttest(dpz_rest);
text(1,0.075,sprintf('p = %.2f',p(1)),'horizontalalignment','center')
text(2,0.075,sprintf('p = %.2f',p(2)),'horizontalalignment','center')
title('Mean FC (Atx vs. Pbo): Rest'); 
set(gca,'xtick',[1 2],'xticklabel',{'Block1';'Block2'})
% --------------------
% ATOMOXETINE RESULTS DURING TASK
% --------------------
subplot(2,2,3); hold on

for isubj =1 : 28
  scatter(pos(isubj),atx_task(isubj,1),'w','markerfacecolor',c(isubj,:))
  scatter(1+pos(isubj),atx_task(isubj,2),'w','markerfacecolor',c(isubj,:))
  line([pos(isubj) 1+pos(isubj)],[atx_task(isubj,1) atx_task(isubj,2)],'color',c(isubj,:))
end

line([-0.5 2.5],[0 0],'linestyle',':','color','k')
axis([0.5 2.5 -0.1 0.1]); tp_editplots

[h,p]=ttest(atx_task);
text(1,0.075,sprintf('p = %.2f',p(1)),'horizontalalignment','center')
text(2,0.075,sprintf('p = %.2f',p(2)),'horizontalalignment','center')
title('Mean FC (Atx vs. Pbo): Task'); 
set(gca,'xtick',[1 2],'xticklabel',{'Block1';'Block2'})

% --------------------
% DONEPEZIL RESULTS DURING TASK
% --------------------
subplot(2,2,4); hold on

for isubj =1 : 28
  scatter(pos(isubj),dpz_task(isubj,1),'w','markerfacecolor',c(isubj,:))
  scatter(1+pos(isubj),dpz_task(isubj,2),'w','markerfacecolor',c(isubj,:))
  line([pos(isubj) 1+pos(isubj)],[dpz_task(isubj,1) dpz_task(isubj,2)],'color',c(isubj,:))
end

line([-0.5 2.5],[0 0],'linestyle',':','color','k')
axis([0.5 2.5 -0.1 0.1]); tp_editplots

[h,p]=ttest(dpz_task);
text(1,0.04,sprintf('p = %.2f',p(1)),'horizontalalignment','center')
text(2,0.04,sprintf('p = %.2f',p(2)),'horizontalalignment','center')
title('Mean FC (Dpz vs. Pbo): Task'); 
set(gca,'xtick',[1 2],'xticklabel',{'Block1';'Block2'})

%% PLOT ALTERED CORRELATIONS

para.alpha = 0.05;
para.nfreq = 1:21;

altered_corr_block1 = pupmod_compute_altered_correlations(fc(:,:,:,:,:,:,1), para)
altered_corr_block2 = pupmod_compute_altered_correlations(fc(:,:,:,:,:,:,2), para)
%%


figure_w;

lims = [0 0.25 0.5];
lims_lab = num2cell([0 25 50]);
lims_ctx = [-0.5 0 0.5; -0.5 0 0.5];
nfreq = 21;

subplot(4,2,3); hold on
plot(altered_corr_block1.n_p_atx(:,2),'r-','linewidth',2)
plot(altered_corr_block1.n_n_atx(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])

subplot(4,2,5); hold on
plot(altered_corr_block1.n_p_atx(:,1),'r-','linewidth',2)
plot(altered_corr_block1.n_n_atx(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])


subplot(4,2,4); hold on
plot(altered_corr_block1.n_p_dpz(:,2),'r-','linewidth',2)
plot(altered_corr_block1.n_n_dpz(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])

subplot(4,2,6); hold on
plot(altered_corr_block1.n_p_dpz(:,1),'r-','linewidth',2)
plot(altered_corr_block1.n_n_dpz(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])

figure_w;

lims = [0 0.25 0.5];
lims_lab = num2cell([0 25 50]);
lims_ctx = [-0.5 0 0.5; -0.5 0 0.5];
nfreq = 21;

subplot(4,2,1); hold on
plot(altered_corr_block2.n_p_atx(:,2),'r-','linewidth',2)
plot(altered_corr_block2.n_n_atx(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])

subplot(4,2,3); hold on
plot(altered_corr_block2.n_p_atx(:,1),'r-','linewidth',2)
plot(altered_corr_block2.n_n_atx(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])


subplot(4,2,2); hold on
plot(altered_corr_block2.n_p_dpz(:,2),'r-','linewidth',2)
plot(altered_corr_block2.n_n_dpz(:,2),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])

subplot(4,2,4); hold on
plot(altered_corr_block2.n_p_dpz(:,1),'r-','linewidth',2)
plot(altered_corr_block2.n_n_dpz(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) lims(1)-0.05 lims(end)])
tp_editplots
pos(3,:)=get(gca,'Position')
axis([0 nfreq lims(1)-0.05 lims(end)])


subplot(4,2,5); hold on
plot(altered_corr_block2.n_p_atx(:,2)-altered_corr_block2.n_p_atx(:,1),'r-','linewidth',2)
plot(altered_corr_block2.n_n_atx(:,2)-altered_corr_block2.n_n_atx(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task-Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) -0.2 0.2])
tp_editplots
pos(3,:)=get(gca,'Position')

subplot(4,2,6); hold on
plot(altered_corr_block2.n_p_dpz(:,2)-altered_corr_block2.n_p_dpz(:,1),'r-','linewidth',2)
plot(altered_corr_block2.n_n_dpz(:,2)-altered_corr_block2.n_n_dpz(:,1),'b-','linewidth',2)
set(gca,'tickdir','out','xtick',[1 5 9 13 17 21],'xticklabel',num2cell([2 4 8 16 32 64]))
set(gca,'tickdir','out','ytick',lims,'yticklabel',lims_lab)
ylabel('Altered corr. [%]')
title('Task-Rest')
axis([0 size(altered_corr_block1.n_p_atx,1) -0.2 0.2])
tp_editplots
pos(3,:)=get(gca,'Position')
% axis([0 nfreq -0.2 0.2])


%% CROSS VALIDATION

% ATOMOXETINE
% -----------
% RESTING STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,2,1,6,1),fc(:,:,:,1,1,6,1),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,2,1,6,2)-fc(:,:,isubj,1,1,6,2);
  other_block_res(isubj,2) = mean(tmp(h==1));
end
% RESTING STATE: Test in block1 with block2 connections
h = ttest(fc(:,:,:,2,1,6,2),fc(:,:,:,1,1,6,2),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,2,1,6,1)-fc(:,:,isubj,1,1,6,1);
  other_block_res(isubj,1) = mean(tmp(h==1));
end
% TASK STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,2,2,6,1),fc(:,:,:,1,2,6,1),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,2,2,6,2)-fc(:,:,isubj,1,2,6,2);
  other_block(isubj,2) = mean(tmp(h==1));
end
% TASK STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,2,2,6,2),fc(:,:,:,1,2,6,2),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,2,2,6,1)-fc(:,:,isubj,1,2,6,1);
  other_block(isubj,1) = mean(tmp(h==1));
end

% DONEPEZIL
% -----------
% RESTING STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,3,1,7,1),fc(:,:,:,1,1,7,1),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,3,1,7,2)-fc(:,:,isubj,1,1,7,2);
  other_block_res(isubj,2) = mean(tmp(h==1));
end
% RESTING STATE: Test in block1 with block2 connections
h = ttest(fc(:,:,:,3,1,7,2),fc(:,:,:,1,1,7,2),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,3,1,6,1)-fc(:,:,isubj,1,1,7,1);
  other_block_res(isubj,1) = mean(tmp(h==1));
end
% TASK STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,3,2,6,1),fc(:,:,:,1,2,7,1),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,3,2,7,2)-fc(:,:,isubj,1,2,7,2);
  other_block(isubj,2) = mean(tmp(h==1));
end
% TASK STATE: Test in block2 with block1 connections
h = ttest(fc(:,:,:,3,2,7,2),fc(:,:,:,1,2,7,2),'dim',3);
for isubj = 1 : 28
  tmp = fc(:,:,isubj,3,2,7,1)-fc(:,:,isubj,1,2,7,1);
  other_block(isubj,1) = mean(tmp(h==1));
end

%% COMPARISON ACROSS BLOCKS (PLACEBO FC)
clear r

for cond = 1:3
  for ifoi = 1 : 25
    for isubj = 1 : 28
      fc1 = fc(:,:,isubj,cond,1,ifoi,1);
      fc2 = fc(:,:,isubj,cond,1,ifoi,2);
      
      [r(isubj,ifoi,cond,1) p(isubj,ifoi,1)] = corr(fc1(mask),fc2(mask));
      
      fc1 = fc(:,:,isubj,cond,2,ifoi,1);
      fc2 = fc(:,:,isubj,cond,2,ifoi,2);
      
      [r(isubj,ifoi,cond,2) p(isubj,ifoi,2)] = corr(fc1(mask),fc2(mask));
    end
  end
end



figure; set(gcf,'color','w');
subplot(2,2,1); hold on

plot(squeeze(mean(r))); tp_editplots
axis([0 14 0 0.6]); 
ylabel('Mean correlation'); 
xlabel('Carrier frequency [Hz]'); 
set(gca,'xtick',[1:2:13],'xticklabel',[2 4 8 16 32 64 128])


%% TASK EFFECTS
clear r
for ifoi = 1 : 13
  for isubj = 1 : 28
    tmp_fc = squeeze(fc(:,:,isubj,1,2,ifoi,:)-fc(:,:,isubj,1,1,ifoi,:));
    block1 = tmp_fc(:,:,1);
    block2 = tmp_fc(:,:,2);
    block1 = block1(mask);
    block2 = block2(mask);
    r(isubj,ifoi) = corr(block1,block2);
  end
end

%% DIFFERENCES FOR INDIV SUBJECTS
load redblue.mat
mask = logical(triu(ones(90,90),1));

for isubj = 1 : 28
  
  
  figure_w;
  subplot(1,4,1); 
  pbo = nanmean(cleandat(:,:,isubj,1,2,1:2),6);
  scale = max([abs(min(pbo(mask))) max(pbo(mask))]);
  imagesc(pbo,[-scale scale]); colormap(redblue); axis square 
  title('Rest')

  subplot(1,4,2); 
  pbo = nanmean(cleandat(:,:,isubj,1,2,1:2),6)-nanmean(cleandat(:,:,isubj,1,1,1:2),6);
  scale = max([abs(min(pbo(mask))) max(pbo(mask))]);
  imagesc(pbo,[-scale scale]); colormap(redblue); axis square 
  title('Task-Rest')
  
  subplot(1,4,3)
  pbo=nanmean(cleandat(:,:,isubj,2,1,1:2),6)-nanmean(cleandat(:,:,isubj,1,1,1:2),6);
  scale = max([abs(min(pbo(mask))) max(pbo(mask))]);
  imagesc(pbo,[-0.1 0.1]); colormap(redblue); axis square 
  title('Atx-Rest (Rest)')
  
  subplot(1,4,4)
  pbo=nanmean(cleandat(:,:,isubj,2,2,1:2),6)-nanmean(cleandat(:,:,isubj,1,2,1:2),6);
  scale = max([abs(min(pbo(mask))) max(pbo(mask))]);
  imagesc(pbo,[-0.1 0.1]); colormap(redblue); axis square 
  title('Atx-Rest (Task)')
  
end

%% GROUP SUBJECTS BASED ON EFFECTS DURING REST AND TASK

for isubj = 1 : 28
  
  eff=nanmean(cleandat(:,:,isubj,2,1,1:2),6)-nanmean(cleandat(:,:,isubj,1,1,1:2),6);
  
  REST = nanmean(eff(mask))>0;

  eff=nanmean(cleandat(:,:,isubj,2,2,1:2),6)-nanmean(cleandat(:,:,isubj,1,2,1:2),6);
  
  TASK = nanmean(eff(mask))>0;
  
  if REST==1 && TASK==1
    group(isubj)=1;
  elseif REST==0 && TASK==1
    group(isubj)=2;
  elseif REST==1 && TASK==0
    group(isubj)=3;
  else REST==0 && TASK==0
    group(isubj)=4;
  end

end
  

eff1 = nanmean(nanmean(cleandat(:,:,group==4,1,1,1:2),6),3)-nanmean(nanmean(cleandat(:,:,group==1,1,1,1:2),6),3);
eff2 = nanmean(nanmean(cleandat(:,:,group==4,1,2,1:2),6),3)-nanmean(nanmean(cleandat(:,:,group==1,1,2,1:2),6),3);


figure_w; 
subplot(1,2,1);
imagesc(eff1,[-0.05 0.05]); colormap(redblue); axis square
subplot(1,2,2);
imagesc(eff2,[-0.05 0.05]); colormap(redblue); axis square

%% POWER (LOAD SIGNAL FROM pupmod_sens_peakfreq.m)
foi_range       = 2.^[1:.25:6]; 

 for ifoi = 1 : 21
   [wavelet,~,opt] = tp_mkwavelet(11.3137,0.5,400);
   
   [~,f,~] = tp_mkwavelet(foi_range(ifoi),0.5,400);
   idx=out.pwelchfreq>f(1)&out.pwelchfreq<f(2);
   fpow(:,ifoi) = squeeze(nanmean(pow(idx,:,1,1),1));
   
   
 end

fc_avg = squeeze(nanmean(nanmean(nanmean(fc(:,:,:,1,1,:,:),1),2),7));
for isubj = 1 : 28
  
  figure_w;
  subplot(1,2,1);
  plot(5:17,fpow(isubj,5:17)); axis square
  set(gca,'XTick',5:4:17,'XTickLabel',num2cell([round(foi_range(5:4:17))]))
set(gca,'ydir','normal'); axis([5 17 min(fpow(:)) max(fpow(:))])
tp_editplots; 
  subplot(1,2,2);
  plot(5:17,fc_avg(isubj,5:17)); axis square
  set(gca,'XTick',5:4:17,'XTickLabel',num2cell([round(foi_range(5:4:17))]))
set(gca,'ydir','normal'); axis([5 17 0 0.2])
tp_editplots
  
end

%%
[r,p]=corr(fc_avg(:,1:20),fpow(:,1:20));

xlabel('Mean Power'); ylabel('Mean FC (OPEC)');
set(gca,'XTick',1:2:20,'XTickLabel',num2cell([round(foi_range(1:2:20))]))
set(gca,'YTick',1:2:20,'YTickLabel',num2cell([round(foi_range(1:2:20))]))
set(gca,'ydir','normal')
colormap(redblue)

%% PERCENT CHANGE

atx_rest = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,2,1,1:3,:),7),1),2),6));
atx_task = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,2,2,1:3,:),7),1),2),6));
pbo_task = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,2,1:3,:),7),1),2),6));
pbo_rest = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,1,1:3,:),7),1),2),6));

dpz_rest = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,3,1,3,:),7),1),2),6));
dpz_task = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,3,2,3,:),7),1),2),6));



mean(100*(atx_task-pbo_task)./pbo_task)
mean(100*(atx_rest-pbo_rest)./pbo_rest)


% ans =$


%% PEAK FREQ PHARMA EFFECTS

[~,max_f]=max(pow);
% pf=pf(max_f)

max_f = squeeze(max_f(:,:,1,:));
max_f=freqs(max_f)

atx_fc = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,2,1,1:3,:),1),2),7),6))- squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,1,1:3,:),1),2),7),6));