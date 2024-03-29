%% pupmod_plot_avgfc
% Plot connectivity across freq
% Last changed: 13-11-2018

clear

v = 12;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

grid  = select_chans(sa_meg_template.grid_cortex3000,400);
  
vis_l = [-20 -86 16]./10;
vis_r = [16 -80 26]./10;
smc_l = [-42 -26 54]./10;
smc_r = [38 -32 48]./10;
aud_l = [-54 -22 10]./10;
aud_r = [52 -24 12]./10;

[~,vis_l_idx]=min(sum((repmat(vis_l,[400 1])-grid).^2,2));
[~,vis_r_idx]=min(sum((repmat(vis_r,[400 1])-grid).^2,2));
[~,smc_l_idx]=min(sum((repmat(smc_l,[400 1])-grid).^2,2));
[~,smc_r_idx]=min(sum((repmat(smc_r,[400 1])-grid).^2,2));
[~,aud_l_idx]=min(sum((repmat(aud_l,[400 1])-grid).^2,2));
[~,aud_r_idx]=min(sum((repmat(aud_r,[400 1])-grid).^2,2));

% vis_l-grid
cleandat = pupmod_loadpowcorr(v);

%% PLOT CORRELATION DURING REST

isabs = 0;

% smc_fc = (squeeze(nanmean(cleandat(:,smc_r_idx,:,1,1,:),1))+squeeze(nanmean(cleandat(smc_l_idx,:,:,1,1,:),2)))/2;
% smc_m_fc = squeeze(cleandat(smc_l_idx,smc_r_idx,:,1,1,:));
% ttest(smc_m_fc,smc_fc)
% 
% vis_fc = (squeeze(nanmean(cleandat(:,vis_r_idx,:,1,1,:),1))+squeeze(nanmean(cleandat(vis_l_idx,:,:,1,1,:),2)))/2;
% vis_m_fc = squeeze(cleandat(vis_l_idx,vis_r_idx,:,1,1,:));
% ttest(vis_m_fc,vis_fc)


if ~isabs
 	smc_m_fc = squeeze(nanmean(cleandat(smc_l_idx,smc_r_idx,:,1,1,:),3));
  smc_s_fc = squeeze(std(squeeze(abs(cleandat(smc_l_idx,smc_r_idx,:,1,1,:)))))/sqrt(size(cleandat,3));

  
  vis_m_fc = squeeze(nanmean(cleandat(vis_l_idx,vis_r_idx,:,1,1,:),3));
  vis_s_fc = squeeze(std(squeeze(abs(cleandat(vis_l_idx,vis_r_idx,:,1,1,:)))))/sqrt(size(cleandat,3));
  aud_m_fc = squeeze(nanmean(cleandat(aud_l_idx,aud_r_idx,:,1,1,:),3));
  aud_s_fc = squeeze(std(squeeze(abs(cleandat(aud_l_idx,aud_r_idx,:,1,1,:)))))/sqrt(size(cleandat,3));
else
  smc_m_fc = squeeze(nanmean(abs(cleandat(smc_l_idx,smc_r_idx,:,1,1,:),3)));
  smc_s_fc = squeeze(std(squeeze(abs(cleandat(smc_l_idx,smc_r_idx,:,1,1,:)))))/sqrt(size(cleandat,3));
  vis_m_fc = squeeze(nanmean(abs(cleandat(vis_l_idx,vis_r_idx,:,1,1,:),3)));
  vis_s_fc = squeeze(std(squeeze(abs(cleandat(vis_l_idx,vis_r_idx,:,1,1,:)))))/sqrt(size(cleandat,3));
end

figure; set(gcf,'color','w'); hold on
subplot(3,2,1); hold on
c = cbrewer('qual', 'Paired',10);
shadedErrorBar([],squeeze(vis_m_fc),squeeze(vis_s_fc),{'color',c(1,:),'linewidth',2});
shadedErrorBar([],squeeze(smc_m_fc),squeeze(smc_s_fc),{'color',c(9,:),'linewidth',2});
shadedErrorBar([],squeeze(aud_m_fc),squeeze(aud_s_fc),{'color',c(5,:),'linewidth',2});

tp_editplots;
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.04 0.08 0.12 0.16],'yticklabel',num2cell([0 0.04 0.08 0.12 0.16]))
axis([0 14 -0.01 0.16]); ylabel('Correlation'); xlabel('Carrier frequency [Hz]')


%% PLOT ROI CORRELATIONS IN SPACE

if ~exist('sa_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres; clear sa
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

par = squeeze(abs(cleandat(:,vis_r_idx,:,1,1,7)));
comp_par = squeeze(nanmean(abs(cleandat(:,vis_r_idx,:,1,1,7)),1));
[h,p,~,t] = ttest(par',repmat(comp_par,[1 400]),'alpha',0.000000025);
par = nanmean(par,2);

h = (h.*sign(t.tstat))>0;
par = par'.*h;
par(isnan(par))=1;
par=par(:);

cmap = plasma ;
cmap      = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
para      = [];
para.clim = [-0.05 0.05];
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.5;
tp_plot_surface(par,sa_template,para)

%% P VALUES (FROM PERMUTATION TEST)
%% PERMUTATION TESTING

nperm = 10000;

all_idx1 = randi(2,[size(cleandat,3),nperm]);
fc = squeeze(nanmean(nanmean(cleandat,2),1));

clear d_dpz_cnt_perm d_dpz_res_perm d_atx_cnt_perm d_atx_res_perm

for iperm = 1 : nperm
  
fprintf('Perm #%d\n',iperm);

idx1 = all_idx1(:,iperm);
idx2 = 3-idx1;
clear permdat_cnt1 permdat_res1

for i = 1 : length(idx1)
%   i
  permdat_cnt1(i,1,:) = squeeze(fc(i,idx1(i),2,:));
  permdat_cnt1(i,2,:) = squeeze(fc(i,idx2(i),2,:));
  permdat_res1(i,1,:) = squeeze(fc(i,idx1(i),1,:));
  permdat_res1(i,2,:) = squeeze(fc(i,idx2(i),1,:));
end

d_atx_cnt_perm(iperm,:) = squeeze(nanmean(permdat_cnt1(:,2,:),1))-squeeze(nanmean(permdat_cnt1(:,1,:),1));
d_atx_res_perm(iperm,:) = squeeze(nanmean(permdat_res1(:,2,:),1))-squeeze(nanmean(permdat_res1(:,1,:),1));

idx1 = all_idx1(:,iperm); 
idx2 = 3-idx1; 
idx1(idx1==2)=3; idx2(idx2==2)=3;
clear permdat_cnt1 permdat_res1

for i = 1 : length(idx1)
%   i
  permdat_cnt2(i,1,:) = squeeze(fc(i,idx1(i),2,:));
  permdat_cnt2(i,2,:) = squeeze(fc(i,idx2(i),2,:));
  permdat_res2(i,1,:) = squeeze(fc(i,idx1(i),1,:));
  permdat_res2(i,2,:) = squeeze(fc(i,idx2(i),1,:));
end

d_dpz_cnt_perm(iperm,:) = squeeze(nanmean(permdat_cnt2(:,2,:),1))-squeeze(nanmean(permdat_cnt2(:,1,:),1));
d_dpz_res_perm(iperm,:) = squeeze(nanmean(permdat_res2(:,2,:),1))-squeeze(nanmean(permdat_res2(:,1,:),1));

end

dpz_res_maxd = max(abs(d_dpz_res_perm),[],2);
dpz_cnt_maxd = max(abs(d_dpz_cnt_perm),[],2);
atx_res_maxd = max(abs(d_atx_res_perm),[],2);
atx_cnt_maxd = max(abs(d_atx_cnt_perm),[],2);

for ifoi = 1 : 13
  
  d_atx_cnt = squeeze(mean(fc(:,2,2,ifoi)-fc(:,1,2,ifoi)));
  d_atx_res = squeeze(mean(fc(:,2,1,ifoi)-fc(:,1,1,ifoi)));
  d_dpz_cnt = squeeze(mean(fc(:,3,2,ifoi)-fc(:,1,2,ifoi)));
  d_dpz_res = squeeze(mean(fc(:,3,1,ifoi)-fc(:,1,1,ifoi)));
  
	p_atx_res(ifoi) = 1-sum(abs(d_atx_res)>abs(d_dpz_res_perm(:,ifoi)))/nperm;
  p_atx_cnt(ifoi) = 1-sum(abs(d_atx_cnt)>abs(d_atx_cnt_perm(:,ifoi)))/nperm;
	p_dpz_res(ifoi) = 1-sum(abs(d_dpz_res)>abs(d_dpz_res_perm(:,ifoi)))/nperm;
  p_dpz_cnt(ifoi) = 1-sum(abs(d_dpz_cnt)>abs(d_dpz_cnt_perm(:,ifoi)))/nperm;
  
end


%%

alpha1 = 0.025;
alpha2 = 0.01;
alpha3 = 0.001;

isabs = 0;

figure; set(gcf,'color','w');
c = cbrewer('qual', 'Paired',10);

if ~isabs
  m_fc = squeeze(nanmean(nanmean(nanmean(cleandat),2),3));
  s_fc = squeeze(std(squeeze(nanmean(nanmean(cleandat),2))))/sqrt(size(cleandat,3));
  fc = squeeze(nanmean(nanmean(cleandat(:,:,:,[1:3],[1:2],1:13),2),1));
  m_fc_subj = squeeze(nanmean(nanmean(cleandat),2));
else
  m_fc = squeeze(nanmean(nanmean(nanmean(abs(cleandat)),2),3));
  s_fc = squeeze(std(squeeze(nanmean(nanmean(abs(cleandat)),2))))/sqrt(size(cleandat,3));
  fc = squeeze(nanmean(nanmean(abs(cleandat(:,:,:,[1:3],[1:2],1:13)),2),1)); 
  m_fc_subj = squeeze(nanmean(nanmean(abs(cleandat)),2));
end

subplot(4,2,1); hold on
title('Rest')
% plot(squeeze(m_fc(1,1,:)),'color',c(1,:),'linewidth',2)
% plot(squeeze(m_fc(2,1,:)),'color',c(2,:),'linewidth',2)
shadedErrorBar([],squeeze(m_fc(1,1,:)),squeeze(s_fc(1,1,:)),{'color',[0.8 0.8 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(2,1,:)),squeeze(s_fc(2,1,:)),{'color',c(8,:),'linewidth',2});
% [h1,p1] = ttest(squeeze(fc(:,2,1,:)),squeeze(fc(:,1,1,:))); 
plot(find(p_atx_res<alpha1),repmat(0.06,[1 sum(p_atx_res<alpha1)]),'ko','markersize',7,'markerfacecolor','k')
plot(find(p_atx_res<alpha2),repmat(0.06,[1 sum(p_atx_res<alpha2)]),'ko','markersize',7,'markerfacecolor','w')
plot(find(p_atx_res<alpha3),repmat(0.06,[1 sum(p_atx_res<alpha3)]),'ro','markersize',7,'markerfacecolor','m')
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.02 0.04 0.06],'yticklabel',num2cell([0 0.02 0.04 0.06]))
axis([0 14 -0.01 0.06]); 
ylabel('Amplitude corr.')
tp_editplots
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,3); hold on
title('Task')
shadedErrorBar([],squeeze(m_fc(1,2,:)),squeeze(s_fc(1,2,:)),{'color',[0.8 0.8 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(2,2,:)),squeeze(s_fc(2,2,:)),{'color',c(8,:),'linewidth',2});
[h2,p2] = ttest(squeeze(fc(:,2,2,:)),squeeze(fc(:,1,2,:))); 
plot(find(p_atx_cnt<alpha1),repmat(0.06,[1 sum(p_atx_cnt<alpha1)]),'ko','markersize',7,'markerfacecolor','k')
plot(find(p_atx_cnt<alpha2),repmat(0.06,[1 sum(p_atx_cnt<alpha2)]),'ko','markersize',7,'markerfacecolor','w')
plot(find(p_atx_cnt<alpha3),repmat(0.06,[1 sum(p_atx_cnt<alpha3)]),'ro','markersize',7,'markerfacecolor','m')
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.02 0.04 0.06],'yticklabel',num2cell([0 0.02 0.04 0.06]))
axis([0 14 -0.01 0.06]); tp_editplots
ylabel('Amplitude corr.'); 
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,2); hold on
shadedErrorBar([],squeeze(m_fc(1,1,:)),squeeze(s_fc(1,1,:)),{'color',[0.5 0.5 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(3,1,:)),squeeze(s_fc(3,1,:)),{'color',c(10,:),'linewidth',2});
[h3,p3] = ttest(squeeze(fc(:,3,1,:)),squeeze(fc(:,1,1,:))); 
plot(find(p_dpz_res<alpha1),repmat(0.06,[1 sum(p_dpz_res<alpha1)]),'ko','markersize',7,'markerfacecolor','k')
plot(find(p_dpz_res<alpha2),repmat(0.06,[1 sum(p_dpz_res<alpha2)]),'ko','markersize',7,'markerfacecolor','w')
plot(find(p_dpz_res<alpha3),repmat(0.06,[1 sum(p_dpz_res<alpha3)]),'ro','markersize',7,'markerfacecolor','m')
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.02 0.04 0.06],'yticklabel',num2cell([0 0.02 0.04 0.06]))
axis([0 14 -0.01 0.06]); tp_editplots
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,4); hold on
shadedErrorBar([],squeeze(m_fc(1,2,:)),squeeze(s_fc(1,2,:)),{'color',[0.8 0.8 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(3,2,:)),squeeze(s_fc(3,2,:)),{'color',c(10,:),'linewidth',2});
[h4,p4] = ttest(squeeze(m_fc_subj(:,3,2,:)),squeeze(m_fc_subj(:,1,2,:)));
plot(find(p_dpz_cnt<alpha1),repmat(0.06,[1 sum(p_dpz_cnt<alpha1)]),'ko','markersize',7,'markerfacecolor','k')
plot(find(p_dpz_cnt<alpha2),repmat(0.06,[1 sum(p_dpz_cnt<alpha2)]),'ko','markersize',7,'markerfacecolor','w')
plot(find(p_dpz_cnt<alpha3),repmat(0.06,[1 sum(p_dpz_cnt<alpha3)]),'ro','markersize',7,'markerfacecolor','m')
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[0 0.02 0.04 0.06],'yticklabel',num2cell([0 0.02 0.04 0.06]))
axis([0 14 -0.01 0.06]); tp_editplots
line([0 14],[0 0],'linestyle',':','color',[0.5 0.5 0.5])
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,5); hold on
par1 = squeeze(std((m_fc_subj(:,2,1,:)-m_fc_subj(:,1,1,:))))/sqrt(size(cleandat,3));
par2 = squeeze(std((m_fc_subj(:,2,2,:)-m_fc_subj(:,1,2,:))))/sqrt(size(cleandat,3));
shadedErrorBar([],squeeze(m_fc(2,1,:))-squeeze(m_fc(1,1,:)),par1,{'color',[0.8 0.8 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(2,2,:))-squeeze(m_fc(1,2,:)),par2,{'color',c(8,:),'linewidth',2});
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[-0.01 0 0.01],'yticklabel',num2cell([-0.01 0 0.01]))
axis([0 14 -0.015 0.015]); tp_editplots
line([0 14],[0 0],'linestyle',':','color',[0.5 0.5 0.5])
ylabel('Percent change [%]'); 
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,6); hold on
par1 = squeeze(std((m_fc_subj(:,3,1,:)-m_fc_subj(:,1,1,:))))/sqrt(size(cleandat,3));
par2 = squeeze(std((m_fc_subj(:,3,2,:)-m_fc_subj(:,1,2,:))))/sqrt(size(cleandat,3));
shadedErrorBar([],squeeze(m_fc(3,1,:))-squeeze(m_fc(1,1,:)),par1,{'color',[0.8 0.8 0.8],'linewidth',2});
shadedErrorBar([],squeeze(m_fc(3,2,:))-squeeze(m_fc(1,2,:)),par2,{'color',c(10,:),'linewidth',2});
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[-0.01 0 0.01],'yticklabel',num2cell([-0.01 0 0.01]))
axis([0 14 -0.015 0.015]); tp_editplots
line([0 14],[0 0],'linestyle',':','color',[0.5 0.5 0.5])
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,7); hold on
par1 = squeeze(std((m_fc_subj(:,2,1,:)-m_fc_subj(:,1,1,:))-(m_fc_subj(:,2,2,:)-m_fc_subj(:,1,2,:))))/sqrt(size(cleandat,3));
shadedErrorBar([],(squeeze(m_fc(2,1,:))-squeeze(m_fc(1,1,:)))-(squeeze(m_fc(2,2,:))-squeeze(m_fc(1,2,:))),par1,{'color',[0.8 0.8 0.8],'linewidth',2})
[~,p_ctx_atx]=ttest(squeeze(fc(:,2,1,:))-squeeze(fc(:,1,1,:)),squeeze(fc(:,2,2,:))-squeeze(fc(:,1,2,:)))
plot(find(p_ctx_atx<alpha1),repmat(0.01,[1 sum(p_ctx_atx<alpha1)]),'ko','markersize',7,'markerfacecolor','k')
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[-0.015 0 0.015],'yticklabel',num2cell([-0.015 0 0.015]))
axis([0 14 -0.015 0.015]); tp_editplots
line([0 14],[0 0],'linestyle',':','color',[0.5 0.5 0.5])
xlabel('Carrier frequency [Hz]')
ylabel('Double diff.'); 
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

subplot(4,2,8); hold on
par2 = squeeze(std((m_fc_subj(:,3,1,:)-m_fc_subj(:,1,1,:))-(m_fc_subj(:,3,2,:)-m_fc_subj(:,1,2,:))))/sqrt(size(cleandat,3));
shadedErrorBar([],(squeeze(m_fc(3,1,:))-squeeze(m_fc(1,1,:)))-(squeeze(m_fc(3,2,:))-squeeze(m_fc(1,2,:))),par1,{'color',[0.8 0.8 0.8],'linewidth',2})
[~,p_ctx_dpz]=ttest(squeeze(fc(:,3,1,:))-squeeze(fc(:,1,1,:)),squeeze(fc(:,3,2,:))-squeeze(fc(:,1,2,:)))
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',num2cell([2 4 8 16 32 64 128]))
set(gca,'tickdir','out','ytick',[-0.015 0 0.015],'yticklabel',num2cell([-0.015 0 0.015]))
axis([0 14 -0.015 0.015]); tp_editplots
line([0 14],[0 0],'linestyle',':','color',[0.5 0.5 0.5])
xlabel('Carrier frequency [Hz]')
a=get(gca,'Position'); a(4) = 0.0960; set(gca,'Position',a)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_avgfc_v%d.pdf',v))



%% PLOT AVG CORRELATION CHANGES IN AAL RSN

k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

load ~/moddules.mat
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',1))
bcn_lab = aal_labels_bcn;

k = zeros(6,1);
rsn = zeros(6,6,28,3,2,13);

for i = 1 : 6
  for j = 1 : 6
  
  idx1 = find(modz == i);
  idx2 = find(modz == j);
  
  for ii = 1 : length(idx1)
    
%     strcmp(lab(idx1{ii},bcn_lab)
  
    rsn(i,j,:,:,:,:) = squeeze(nanmean(nanmean(cleandat(idx1,idx2,:,:,:,:),1),2));
  
  end
end
end
  %%
figure; set(gcf,'color','w'); hold on

% ATX vs PLACEBO REST
subplot(3,2,1); hold on

m_atx = squeeze(nanmean(nanmean(rsn(:,:,:,2,1,6),3),2));
m_pbo = squeeze(nanmean(nanmean(rsn(:,:,:,1,1,6),3),2));

s_atx = squeeze(nanstd(nanmean(rsn(:,:,:,2,1,6),2),[],3))/sqrt(28);
s_pbo = squeeze(nanstd(nanmean(rsn(:,:,:,1,1,6),2),[],3))/sqrt(28);

bar([1:3:16],m_atx,0.3,'facecolor',[1 0.2 0],'edgecolor',[1 1 1])
bar([2:3:17],m_pbo,0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[1 1 1])
 
idx1 = 1:3:16; idx2 = 2:3:17;
for i = 1 : 6
  line([idx1(i) idx1(i)],[m_atx(i)-s_atx(i) m_atx(i)+s_atx(i)],'color',[0.7 0 0])
  line([idx2(i) idx2(i)],[m_pbo(i)-s_pbo(i) m_pbo(i)+s_pbo(i)],'color',[0.3 0.3 0.3])
end

axis([0 17.5 0 0.075])
tp_editplots; ylabel('Corr. coeff.');
set(gca,'xtick',[1.5:3:16.5],'xticklabel',{'DMN';'SMN';'VIS';'SUB';'IFN';'FPN'})

% ATX vs PLACEBO TASK
subplot(3,2,3); hold on

m_atx = squeeze(nanmean(nanmean(rsn(:,:,:,2,2,6),3),2));
m_pbo = squeeze(nanmean(nanmean(rsn(:,:,:,1,2,6),3),2));

s_atx = squeeze(nanstd(nanmean(rsn(:,:,:,2,2,6),2),[],3))/sqrt(28);
s_pbo = squeeze(nanstd(nanmean(rsn(:,:,:,1,2,6),2),[],3))/sqrt(28);

bar([1:3:16],m_atx,0.3,'facecolor',[1 0.2 0],'edgecolor',[1 1 1])
bar([2:3:17],m_pbo,0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[1 1 1])

idx1 = 1:3:16; idx2 = 2:3:17;
for i = 1 : 6
  line([idx1(i) idx1(i)],[m_atx(i)-s_atx(i) m_atx(i)+s_atx(i)],'color',[0.7 0 0])
  line([idx2(i) idx2(i)],[m_pbo(i)-s_pbo(i) m_pbo(i)+s_pbo(i)],'color',[0.3 0.3 0.3])
end

axis([0 17.5 0 0.075])
tp_editplots; ylabel('Corr. coeff.');
set(gca,'xtick',[1.5:3:16.5],'xticklabel',{'DMN';'SMN';'VIS';'SUB';'IFN';'FPN'})

% ATX vs PLACEBO TASK
subplot(3,2,2); hold on

m_atx = squeeze(nanmean(nanmean(rsn(:,:,:,3,1,7),3),2));
m_pbo = squeeze(nanmean(nanmean(rsn(:,:,:,1,1,7),3),2));

s_atx = squeeze(nanstd(nanmean(rsn(:,:,:,3,1,7),2),[],3))/sqrt(28);
s_pbo = squeeze(nanstd(nanmean(rsn(:,:,:,1,1,7),2),[],3))/sqrt(28);

bar([1:3:16],m_atx,0.3,'facecolor',[0 0.2 1],'edgecolor',[1 1 1])
bar([2:3:17],m_pbo,0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[1 1 1])

idx1 = 1:3:16; idx2 = 2:3:17;
for i = 1 : 6
  line([idx1(i) idx1(i)],[m_atx(i)-s_atx(i) m_atx(i)+s_atx(i)],'color',[0 0 0.7])
  line([idx2(i) idx2(i)],[m_pbo(i)-s_pbo(i) m_pbo(i)+s_pbo(i)],'color',[0.3 0.3 0.3])
end

axis([0 17.5 0 0.05])
tp_editplots; ylabel('Corr. coeff.');
set(gca,'xtick',[1.5:3:16.5],'xticklabel',{'DMN';'SMN';'VIS';'SUB';'IFN';'FPN'})

% ATX vs PLACEBO TASK
subplot(3,2,4); hold on

m_atx = squeeze(nanmean(nanmean(rsn(:,:,:,3,2,7),3),2));
m_pbo = squeeze(nanmean(nanmean(rsn(:,:,:,1,2,7),3),2));

s_atx = squeeze(nanstd(nanmean(rsn(:,:,:,3,2,7),2),[],3))/sqrt(28);
s_pbo = squeeze(nanstd(nanmean(rsn(:,:,:,1,2,7),2),[],3))/sqrt(28);

bar([1:3:16],m_atx,0.3,'facecolor',[0 0.2 1],'edgecolor',[1 1 1])
bar([2:3:17],m_pbo,0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[1 1 1])

idx1 = 1:3:16; idx2 = 2:3:17;
for i = 1 : 6
  line([idx1(i) idx1(i)],[m_atx(i)-s_atx(i) m_atx(i)+s_atx(i)],'color',[0 0 0.7])
  line([idx2(i) idx2(i)],[m_pbo(i)-s_pbo(i) m_pbo(i)+s_pbo(i)],'color',[0.3 0.3 0.3])
end

axis([0 17.5 0 0.05])
tp_editplots; ylabel('Corr. coeff.');
set(gca,'xtick',[1.5:3:16.5],'xticklabel',{'DMN';'SMN';'VIS';'SUB';'IFN';'FPN'})

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_RSN_bar.pdf'))

%% PLOT CHANGES IN CONNECTIVITY ON SURFACE

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

close all

ifoi = 7; icond = 1;

cmap = plasma;

fc = squeeze(nanmean(cleandat(:,:,:,:,:,ifoi),2));
[~,p,~,t] = ttest(fc(:,:,3,icond),fc(:,:,1,icond),'dim',2);
par = t.tstat;
% par = mean((fc(:,:,3,icond)-fc(:,:,1,icond))./fc(:,:,1,icond),2);
% par(outp_atx.p_cnt_p_atx(:,icond,ifoi)>=0.025) = 0;
load redblue.mat
cmap = redblue;
para = [];
para.clim = [-5 5];
para.cmap = cmap;
para.grid = grid;
para.dd = .75;
tp_plot_surface(par,sa_template,para)

% print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_atx_f%d_c%d_v%d.pdf',ifoi,icond,v))

% ifoi = 7; icond = 1;
% 
% cmap = viridis;
% 
% par = emp.n_n_dpz_pervoxel(:,ifoi,icond);
% par(outp_dpz.n_cnt_p_dpz(:,icond,ifoi)>=0.025) = 0;
% 
% cmap = [cmap(50:end-15,:); 0.98*ones(1,3); cmap(50:end-15,:)];
% para = [];
% para.clim = [-0.65 0.65];
% para.cmap = cmap;
% para.grid = grid;
% para.dd = 0.75;
% tp_plot_surface(par,sa_template,para)

%% CORRELATIONS PER SUBJECT
for ifoi = 1 : 13
  fc = squeeze(nanmean(cleandat(:,:,:,:,:,ifoi),2));
  for isubj = 1 : 28
    larger_than_pbo(isubj,2,ifoi) = 100*sum(fc(:,isubj,2,2)>fc(:,isubj,1,2))/400;
    larger_than_pbo(isubj,1,ifoi) = 100*sum(fc(:,isubj,2,1)>fc(:,isubj,1,1))/400;
    r_atx(isubj,1) = corr(fc(:,isubj,2,1),fc(:,isubj,1,1));
    r_atx(isubj,2) = corr(fc(:,isubj,2,2),fc(:,isubj,1,2));
  end
end

%%
figure; set(gcf,'color','w'); hold on
m = squeeze(mean(larger_than_pbo));
s = squeeze(std(larger_than_pbo,[],1))/sqrt(28);
shadedErrorBar(1:13,m(1,:),s(1,:),{'color',[0.7 0.2 0]})
shadedErrorBar(1:13,m(2,:),s(2,:),{'color',[0.3 0 0]})
line([0 14],[50 50],'linestyle',':')
axis([0 14 0 100]); tp_editplots;
xlabel('Frequency [Hz]'); ylabel('Fraction of voxels larger than Pbo.')
%%
fc = squeeze(nanmean(cleandat(:,:,:,:,:,7),2));
for isubj = 1 : 28
  smaller_than_pbo(isubj,2) = 100*sum(fc(:,isubj,3,2)<fc(:,isubj,1,2))/400;
  smaller_than_pbo(isubj,1) = 100*sum(fc(:,isubj,3,1)<fc(:,isubj,1,1))/400;
end

%% PLOT ALL FC CONTRASTS
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip');
cmap = cmap(end:-1:1,:);
load redblue.mat

for ifoi = 1 : 13
  
  figure; set(gcf,'color','w');
  
  subplot(2,2,1); hold on
  title('Rest - ATX')

  [h,~,~,t]=ttest(cleandat(:,:,:,2,1,ifoi),cleandat(:,:,:,1,1,ifoi),'dim',3);
  
  imagesc(t.tstat.*h,[-1.96 1.96]); colormap(redblue)
  axis equal off
  subplot(2,2,2); hold on
  title('Task - ATX')

  [h,~,~,t]=ttest(cleandat(:,:,:,2,2,ifoi),cleandat(:,:,:,1,2,ifoi),'dim',3);
  
  imagesc(t.tstat.*h,[-1.96 1.96]); colormap(redblue)
  axis equal off
  
  subplot(2,2,3); hold on
  title('Rest - DPZ')

  [h,~,~,t]=ttest(cleandat(:,:,:,3,1,ifoi),cleandat(:,:,:,1,1,ifoi),'dim',3);
  
  imagesc(t.tstat.*h,[-1.96 1.96]); colormap(redblue)
  axis equal off
  
  subplot(2,2,4); hold on
  title('Task - DPZ')
  
  [h,~,~,t]=ttest(cleandat(:,:,:,3,2,ifoi),cleandat(:,:,:,1,2,ifoi),'dim',3);
  
  imagesc(t.tstat.*h,[-1.96 1.96]); colormap(redblue)
  axis equal off
  
  set(gcf,'Renderer','Painters')
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_fcmatrices_f%d.pdf',ifoi))
  
end



