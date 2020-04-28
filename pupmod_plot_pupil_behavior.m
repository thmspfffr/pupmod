%% PLOT PUPIL AND BEHAVIOR

% Content:

% 1. CORTEX
%---------------
%   1.1. Fraction of altered correlations
%   1.2. Drug effect on FC
%   1.3. Correlations of behavior with FC (per edge)
%   1.4. Task vs. Rest
%---------------
% 2. VTPM (40x40)
%   1.1. Fraction of altered correlations
%   1.2. Drug effect on FC
%   1.3. Correlations of behavior with FC (per edge)
%   1.4. Task vs. Rest
%---------------

%%
clear

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pconn/matlab
addpath ~/pupmod/matlab

% v  = 23;
% cleandat = pupmod_loadpowcorr(v,SUBJLIST,0);
% fc = cat(4,fc(:,:,:,:,:,:,1),fc(:,:,:,:,:,:,2));
% fc = (fc-nanmean(fc,4))./nanstd(fc,[],4);
% fc = cat(7,fc(:,:,:,1:3,:,:),fc(:,:,:,4:6,:,:));
%
para.str_behav  = 'count';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt       = behav;

% tmp = reshape(behav_cnt,[3*2 28]);
% tmp = (tmp-nanmean(tmp,1))./nanstd(tmp,1);
% behav_cnt = reshape(tmp,[3 28 2]);

para.str_behav  = 'numb_switches';
behav           = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn      = behav;
behav_bttn      = permute(behav_bttn,[2 1 3]);
%
behav_pooled = (behav_cnt+behav_bttn)./2;
behav_pooled(isnan(behav_bttn))=behav_cnt(isnan(behav_bttn));
behav_pooled(isnan(behav_cnt))=behav_bttn(isnan(behav_cnt));

para= [];
para.t = [3000 7000];
pup = pp_loadpupil(para);
% pup(:,:,:,1) = rest | :,2) = counting | :,3) = butt

pup_pooled = nanmean(pup(:,:,:,1:2),4);
% pup_pooled(20,:,:)= [];

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%%
cols = cbrewer('seq', 'Reds', 28);

figure; set(gcf,'color','w');
subplot(3,3,1:2); hold on

col_rand = (rand-0.5)/2;

par = squeeze(nanmean(pup,3));
par_pooled = squeeze(nanmean(pup_pooled,3));

r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  
  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
%   plot(r(isubj)+6,par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
%   plot(r(isubj)+8,par(isubj,2,2),'.','markersize',10,'color',cols(isubj,:))
%   line([r(isubj)+6 r(isubj)+8],[par(isubj,1,2) par(isubj,2,2)],'color',cols(isubj,:))
end

plot(0,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(2,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(4,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')


axis([-2 25 4000 12000]); tp_editplots

subplot(3,3,4:5); hold on

% col_rand = (rand-0.5)/2;

par_pooled = squeeze(nanmean(behav_pooled,3))';

% r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  
  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
%   plot(r(isubj)+6,par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
%   plot(r(isubj)+8,par(isubj,2,2),'.','markersize',10,'color',cols(isubj,:))
%   line([r(isubj)+6 r(isubj)+8],[par(isubj,1,2) par(isubj,2,2)],'color',cols(isubj,:))
end

plot(0,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(2,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(4,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')


axis([-02 25 0 150]); tp_editplots; %lsline


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_pupil_behavior.pdf'))

%% SEPARATELY FOR REST AND TASK BLOCKS


cols = cbrewer('seq', 'Reds', 28);

figure; set(gcf,'color','w');
subplot(3,3,1:2); hold on

col_rand = (rand-0.5)/2;

% ------------------------------------
% REST PUPIL FIRST
% ------------------------------------
pup_rest = nanmean(pup(:,:,:,1),3);
par_pooled = squeeze(nanmean(pup_rest,3));

r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
end

plot(0,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(2,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(4,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

axis([-2 25 4000 12000]); tp_editplots

% ------------------------------------
% COOUNTING PUPIL FIRST
% ------------------------------------

pup_cnt = nanmean(pup(:,:,:,2),3);
par_pooled = squeeze(nanmean(pup_cnt,3));

r=rand(28,1)-0.5;

for isubj=1:28
  col_rand = (rand-0.5)/5;

  plot(r(isubj)+6,par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+8,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+10,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  line([r(isubj)+6 r(isubj)+8],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+8 r(isubj)+10],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
end

plot(6,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(8,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(10,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([6 8],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([8 10],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')

subplot(3,3,4:5); hold on

% col_rand = (rand-0.5)/2;

par_pooled = squeeze(nanmean(permute(behav_cnt,[2 1 3]),3));
par_pooled = par_pooled(~any(isnan(par_pooled),2),:);
% r=rand(28,1)-0.5;

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj),par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+2,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+4,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  
  line([r(isubj) r(isubj)+2],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+2 r(isubj)+4],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
%   plot(r(isubj)+6,par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
%   plot(r(isubj)+8,par(isubj,2,2),'.','markersize',10,'color',cols(isubj,:))
%   line([r(isubj)+6 r(isubj)+8],[par(isubj,1,2) par(isubj,2,2)],'color',cols(isubj,:))
end

plot(0,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(2,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(4,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([0 2],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([2 4],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')


par_pooled = squeeze(nanmean(permute(behav,[2 1 3]),3))';
par_pooled = par_pooled(~any(isnan(par_pooled),2),:);
% r=rand(28,1)-0.5;

for isubj=1:size(par_pooled,1)
  col_rand = (rand-0.5)/5;

  plot(r(isubj)+6,par_pooled(isubj,2,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+8,par_pooled(isubj,1,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  plot(r(isubj)+10,par_pooled(isubj,3,1),'.','markersize',10,'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])

  
  line([r(isubj)+6 r(isubj)+8],[par_pooled(isubj,2,1) par_pooled(isubj,1,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
  line([r(isubj)+8 r(isubj)+10],[par_pooled(isubj,1,1) par_pooled(isubj,3,1)],'color',[0.7+col_rand 0.7+col_rand 0.7+col_rand])
 
%   plot(r(isubj)+6,par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
%   plot(r(isubj)+8,par(isubj,2,2),'.','markersize',10,'color',cols(isubj,:))
%   line([r(isubj)+6 r(isubj)+8],[par(isubj,1,2) par(isubj,2,2)],'color',cols(isubj,:))
end

plot(6,nanmean(par_pooled(:,2)),'.','markersize',20,'color','k')
plot(8,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(10,nanmean(par_pooled(:,3)),'.','markersize',20,'color','k')

line([6 8],[nanmean(par_pooled(:,2)) nanmean(par_pooled(:,1))],'color','k')
line([ 10],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,3))],'color','k')



axis([-02 25 0 150]); tp_editplots; %lsline


print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_pupil_behavior_resttaskseparately.pdf'))


%%

cond = [2 3];
r=rand(28,1)-0.5;
par = squeeze(nanmean(pup,3));
par([14 20],:,:)=[];

for icond = 1 : 2

figure; set(gcf,'color','w');
subplot(3,3,1:2); hold on




c = cond(icond);
if icond == 1
  cols = cbrewer('seq', 'Reds', 28);

else
  cols = cbrewer('seq', 'Blues', 28);

end

for isubj=1:size(par,1)
  plot(r(isubj),par(isubj,1,1),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+2,par(isubj,c,1),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj) r(isubj)+2],[par(isubj,1,1) par(isubj,c,1)],'color',cols(isubj,:))
  
  plot(r(isubj)+6,par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+8,par(isubj,c,2),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj)+6 r(isubj)+8],[par(isubj,1,2) par(isubj,c,2)],'color',cols(isubj,:))
end
for isubj=1:size(par,1)
  if any(isnan(par_pooled(isubj,:))); continue; end
  plot(r(isubj)+12,par_pooled(isubj,1),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+14,par_pooled(isubj,c),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj)+12 r(isubj)+14],[par_pooled(isubj,1) par_pooled(isubj,c)],'color',cols(isubj,:))
end


plot(-1,nanmean(par(:,1,1)),'.','markersize',20,'color','k')
plot(3,nanmean(par(:,c,1)),'.','markersize',20,'color','k')
line([-1 3],[nanmean(par(:,1,1)) nanmean(par(:,c,1))],'color','k')

plot(5,nanmean(par(:,1,2)),'.','markersize',20,'color','k')
plot(9,nanmean(par(:,c,2)),'.','markersize',20,'color','k')
line([5 9],[nanmean(par(:,1,2)) nanmean(par(:,c,2))],'color','k')

plot(11,nanmean(par_pooled(:,1)),'.','markersize',20,'color','k')
plot(15,nanmean(par_pooled(:,c)),'.','markersize',20,'color','k')
line([11 15],[nanmean(par_pooled(:,1)) nanmean(par_pooled(:,c))],'color','k')

axis([-2 25 4000 12000]); tp_editplots
ylabel('Pupil diameter [a.U.]');

subplot(3,3,3); hold on
for isubj=1:size(par,1)
plot(par(isubj,3,1)-par(isubj,1,1),par(isubj,c,2)-par(isubj,1,2),'.','markersize',10,'color',cols(isubj,:))
end

plot(nanmean(par(:,3,1)-par(:,1,1)),nanmean(par(:,c,2)-par(:,1,2)),'.','markersize',20,'color','k')


[rr,p]=corr(par(:,3,1)-par(:,1,1),par(:,c,2)-par(:,1,2));
text(-3000, 2750,sprintf('r=%3f\np=%3f',rr,p),'fontsize',6)

axis([-3500 3500 -3500 3500]); tp_editplots; %lsline
axis square
line([-3000 3000],[-3000 3000],'color','k'); %lsline
line([0 0],[-3000 3000],'color',[0.5 0.5 0.5],'linestyle',':'); %lsline
line([-3000 3000],[0 0],'color',[0.5 0.5 0.5],'linestyle',':'); 
xlabel('\DeltaDiameter (Rest)');
ylabel('\DeltaDiameter (Task)');
set(gca,'tickdir','out','xtick',[-3000 0 3000],'xticklabel',num2str([-3000; 0; 3000]));
set(gca,'tickdir','out','ytick',[-3000 0 3000],'yticklabel',num2str([-3000; 0; 3000]));
% set(gca,'XTick'

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmpd_pupil_cond%d.eps',icond))

end



%%
subplot(3,3,4:5); hold on

par=nanmean(behav_cnt,3)';

for isubj=1:size(par,1)
  plot(r(isubj),par(isubj,1),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+2,par(isubj,3),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj) r(isubj)+2],[par(isubj,1,1) par(isubj,3,1)],'color',cols(isubj,:))

end


plot(-1,mean(par(:,1)),'.','markersize',20,'color','k')
plot(3,mean(par(:,3)),'.','markersize',20,'color','k')
line([-1 3],[mean(par(:,1)) mean(par(:,3))],'color','k')

par=nanmean(behav_bttn,3)';

for isubj=1:28
  plot(r(isubj)+6,par(isubj,1),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+8,par(isubj,3),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj)+6 r(isubj)+8],[par(isubj,1,1) par(isubj,3,1)],'color',cols(isubj,:))

end

plot(5,nanmean(par(:,1)),'.','markersize',20,'color','k')
plot(9,nanmean(par(:,3)),'.','markersize',20,'color','k')
line([5 9],[nanmean(par(:,1)) nanmean(par(:,3))],'color','k')

par=nanmean(behav_pooled,3)';

for isubj=1:28
  plot(r(isubj)+12,par(isubj,1),'.','markersize',10,'color',cols(isubj,:))
  plot(r(isubj)+14,par(isubj,3),'.','markersize',10,'color',cols(isubj,:))
  line([r(isubj)+12 r(isubj)+14],[par(isubj,1) par(isubj,3)],'color',cols(isubj,:))

end

plot(11,nanmean(par(:,1)),'.','markersize',20,'color','k')
plot(15,nanmean(par(:,3)),'.','markersize',20,'color','k')
line([11 15],[nanmean(par(:,1)) nanmean(par(:,3))],'color','k')

axis([-2 25 00 160]); tp_editplots
ylabel('Percieved switches');
set(gca,'tickdir','out','ytick',[0 20 40 60 80 100 120 140 160],'xticklabel',num2str([0; 20; 40 ;60 ;80; 100; 120; 140; 160]));


clear par
par1 = nanmean(behav_cnt,3)';
par2 = nanmean(behav_bttn,3)';

[i2,j2]=find(isnan(par2));
par1(i2,:)=[]; par2(i2,:)=[];

subplot(3,3,6); hold on
for isubj=1:size(par1,1)
plot(par1(isubj,3)-par1(isubj,1),par2(isubj,3)-par2(isubj,1),'.','markersize',10,'color',cols(isubj,:))
end

[rr,p]=corr(par1(:,3)-par1(:,1),par2(:,3)-par2(:,1));
text(-50, 50,sprintf('r=%3f\np=%3f',rr,p),'fontsize',6)

plot(mean(par1(:,3)-par1(:,1)),mean(par2(:,3)-par2(:,1)),'.','markersize',20,'color','k')

axis([-60 60 -60 60]); tp_editplots; %lsline
line([-50 50],[-50 50],'color','k'); %lsline
axis square
line([0 0],[-50 50],'color',[0.5 0.5 0.5],'linestyle',':'); %lsline
line([-50 50],[0 0],'color',[0.5 0.5 0.5],'linestyle',':'); 
xlabel('\DeltaSwitches (Counting)');
ylabel('\DeltaSwitches (Button)');
set(gca,'tickdir','out','xtick',[-50 0 50],'xticklabel',num2str([-50; 0; 50]));
set(gca,'tickdir','out','ytick',[-50 0 50],'yticklabel',num2str([-50; 0; 50]));
% set(gca,'XTick'

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_pupil_behavior_dpz.pdf'))