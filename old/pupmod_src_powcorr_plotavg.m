v = 1;

%%


load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

for icond = 1 : 2

    s_fc(:,:,1,icond) = squeeze(nanmean(nanmean(cleandat(:,:,:,1,icond,:))));
   	s_fc(:,:,2,icond) = squeeze(nanmean(nanmean(cleandat(:,:,:,2,icond,:))));
    s_fc(:,:,3,icond) = squeeze(nanmean(nanmean(cleandat(:,:,:,3,icond,:))));


end
%% STATS

para = [];
para.paired = 1;
para.nperm = 10000;
para.tail = 0;
para.alpha = 0.05;

for icond = 1 : 2
  for ifoi = 1 : 13
    
    [h(1,icond,ifoi),p(1,icond,ifoi)]=tp_permtest(s_fc(:,ifoi,1,icond),s_fc(:,ifoi,2,icond),para);
    [h(2,icond,ifoi),p(2,icond,ifoi)]=tp_permtest(s_fc(:,ifoi,1,icond),s_fc(:,ifoi,3,icond),para);
    [h(3,icond,ifoi),p(3,icond,ifoi)]=tp_permtest(s_fc(:,ifoi,2,icond),s_fc(:,ifoi,3,icond),para);
    
  end
end

  
  

%% PLOT
f = unique(round(2.^[1:.5:7]));

figure; hold on

subplot(2,2,1); hold on

s1 = std(s_fc(:,:,1,1))./sqrt(28);
s2 = std(s_fc(:,:,2,1))./sqrt(28);

% plot(log10(f),squeeze(mean(s_fc(:,:,1,1))),'linewidth',2,'color',[0.7 0.7 0.7]);
% plot(log10(f),squeeze(mean(s_fc(:,:,2,1))),'linewidth',2,'color','r');
shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,1,1))),s1,{'linewidth',2,'color',[0.7 0.7 0.7]});
shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,2,1))),s2,{'linewidth',2,'color',[1 0.3 0.1]});

plot(log10(f(find(squeeze(h(1,1,:))))),ones(sum(squeeze(h(1,1,:))),1).*squeeze(mean(s_fc(:,find(squeeze(h(1,1,:))),1,1))),'.','markersize',20,'color','k');

xlabel('Frequency [Hz]'); ylabel('Amplitude correlation');
title('Rest - Atomoxetine')


axis([-0.05 2.2 -0.005 0.04])
set(gca,'tickdir','out','xtick',log10(f(1:2:end)),'xticklabel',f(1:2:end))
subplot(2,2,3); hold on

s1 = std(s_fc(:,:,1,2))./sqrt(28);
s2 = std(s_fc(:,:,2,2))./sqrt(28);

shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,1,2))),s1,{'linewidth',2,'color',[0.7 0.7 0.7]});
shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,2,2))),s2,{'linewidth',2,'color',[1 0.3 0.1]});

plot(log10(f(find(squeeze(h(1,2,:))))),ones(sum(squeeze(h(1,2,:))),1).*0,'.','markersize',20,'color','k');

axis([-0.05 2.2 -0.005 0.04])
set(gca,'tickdir','out','xtick',log10(f(1:2:end)),'xticklabel',f(1:2:end))

xlabel('Frequency [Hz]'); ylabel('Amplitude correlation')
title('Task - Atomoxetine')

subplot(2,2,2); hold on
set(gca,'tickdir','out','xtick',log10(f(1:2:end)),'xticklabel',f(1:2:end))

s1 = std(s_fc(:,:,1,1))./sqrt(28);
s2 = std(s_fc(:,:,3,1))./sqrt(28);

shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,1,1))),s1,{'linewidth',2,'color',[0.7 0.7 0.7]});
shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,3,1))),s2,{'linewidth',2,'color',[0.2 0.5 0.9]});

plot(log10(f(find(squeeze(h(2,1,:))))),ones(sum(squeeze(h(2,1,:))),1).*0,'.','markersize',20,'color','k');

xlabel('Frequency [Hz]'); ylabel('Average amplitude correlation')
title('Rest - Donepezil')

axis([-0.05 2.2 -0.005 0.04])
set(gca,'tickdir','out','xtick',log10(f(1:2:end)),'xticklabel',f(1:2:end))

subplot(2,2,4); hold on

s1 = std(s_fc(:,:,1,2))./sqrt(28);
s2 = std(s_fc(:,:,3,2))./sqrt(28);

shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,1,2))),s1,{'linewidth',2,'color',[0.7 0.7 0.7]});
shadedErrorBar(log10(f),squeeze(mean(s_fc(:,:,3,2))),s2,{'linewidth',2,'color',[0.2 0.5 0.9]});

plot(log10(f(find(squeeze(h(2,2,:))))),ones(sum(squeeze(h(2,2,:))),1).*0,'.','markersize',20,'color','k');

xlabel('Frequency [Hz]'); ylabel('Average amplitude correlation')
title('Task - Donepezil')

axis([-0.05 2.2 -0.005 0.04])
set(gca,'tickdir','out','xtick',log10(f(1:2:end)),'xticklabel',f(1:2:end))

%% AAL grid
cnt = 0; clear a 
for i = 1 : 3000
 if isempty(aalgrid.labels{i})
  continue
 else
   cnt = cnt + 1;
   a = aalgrid.labels{i};
   
   idx = find(strcmp(aalgrid.labels,a));
   
   for ik = 1 : length(idx)
     
     
    aalgrid.labels{idx(ik)} = [];
    
   end
   
   labels{cnt} = a;
   
   
 end
end
   

%% PLOT ON SURFACE OF BRAIN

load sa_meg_template;

grid  = sa_meg_template.grid_cortex3000;
g1    = sa_meg_template.grid_cortex3000;
g2    = sa_meg_template.cortex10K.vc;
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = sa_meg_template.grid_cortex3000;

load aalmask_grid_cortex3000.mat
% 
% fc_atx_surf(:,1) = nanmean(fc_atx(:,:,1),2);
% fc_atx_surf(:,2) = nanmean(fc_atx(:,:,2),2);
% 
% % fc_dpz_surf(:,1) = nanmean(fc_dpz(:,:,1),2);
% fc_dpz_surf(:,2) = nanmean(fc_dpz(:,:,2),2);
% % 
par2plot_atx = zeros(3000,1);

for i = 3 : 3

  par2plot_atx(find(aalgrid.mask==i),1) =1;
%   par2plot_atx(find(aalgrid.mask==i),2) = fc_atx_surf(i,2);

%   par2plot_dpz(find(aalgrid.mask==i),1) = fc_dpz_surf(i,1);
%   par2plot_dpz(find(aalgrid.mask==i),2) = fc_dpz_surf(i,2);

%   par2plot_atx_dd(find(aalgrid.mask==i),1) = fc_atx_surf(i,2) - fc_atx_surf(i,1);
%   par2plot_dpz_dd(find(aalgrid.mask==i),1) = fc_dpz_surf(i,2) - fc_dpz_surf(i,1);

end


para = [];
para.colorlimits  = [-1 1];
par_interp        = spatfiltergauss(par2plot_atx,g1,dd,g2);

tp_showsource(par_interp,hot,sa_meg_template,para)



