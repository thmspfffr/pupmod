%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

% Last changed: 13-11-2018

clear

v = 20;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

fprintf('Loading grid...\n')
grid  = select_chans(sa_meg_template.grid_cortex3000,400); fprintf('Loading grid... Done\n')

cleandat = pupmod_loadpowcorr(v,0);

if v == 20
  tmp = tp_create_grid('vtpm');
  clim = [-0.025 0.025];
  idx = 1:46;
  idx = ~ismember(idx,[21 22 23 44 45 46]);
  reg = tmp.tissuelabel_4mm(idx);

end
%% PLOT FC MATRIX FOR DIFFERENT CONDITIONS
FOI = 6:7

if v == 12
  for ifoi = FOI

    figure; set(gcf,'color','w')
    cond = [1 2 3];
    [~,front_to_back] = sort(grid(:,2),'descend');
    left  = find(grid(:,1)<0);
    right = find(grid(:,1)>0);
    for i = 1 : 3
      fc = nanmean(nanmean(cleandat(:,:,:,cond(i),1,ifoi,:),7),3);

      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));

      fc_rest = [fc1 fc2; fc3 fc4];

      fc = nanmean(nanmean(cleandat(:,:,:,cond(i),2,ifoi,:),7),3);

      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));

      fc_task = [fc1 fc2; fc3 fc4];

      fc = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(1,3,i); imagesc(fc,[0.02 0.1]); axis square off
      colormap(inferno)

    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_f%d_v%d.pdf',ifoi,v))
  end
elseif v== 20
  for ifoi = FOI

    figure; set(gcf,'color','w')
    cond = [1 2 3];
    
    for i = 1 : 3
      
      fc_rest = nanmean(nanmean(cleandat(idx,idx,:,cond(i),1,ifoi,:),7),3);
      fc_task = nanmean(nanmean(cleandat(idx,idx,:,cond(i),2,ifoi,:),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
%       if ifoi ~= 7
      subplot(1,3,i); imagesc(fc,[0.04 0.12]); axis square %off
      
      set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
      colormap(inferno); tp_editplots; set(gca,'FontSize',5)

    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_f%d_v%d.pdf',ifoi,v))
  end
end
 
%% DRUG CONTRASTS: AVERAGED ACROSS BLOCKS

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
cond = [2 3];

SUBJ = 1:28;
FOI = 6:7;
CLIM = [-0.02 0.02];

if v == 12
  
  [~,front_to_back] = sort(grid(:,2),'descend');
  
  left  = find(grid(:,1)<0);
  right = find(grid(:,1)>0);
  for ifoi = FOI
    
    figure; set(gcf,'color','w')
    
    for i = 1 : 2
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),1,ifoi,:),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,1,ifoi,:),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_rest = [fc1 fc2; fc3 fc4];
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),2,ifoi,:),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,2,ifoi,:),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_task = [fc1 fc2; fc3 fc4];
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      
      subplot(2,2,i); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Blocks 1&2)'); else; title('Dpz - Pbo (Blocks 1&2)'); end
      colormap(cmap)
      tp_colorbar('Difference in correlation'); 
    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_f%d_v%d.pdf',ifoi,v))

  end
else
  for ifoi = FOI
    figure; set(gcf,'color','w')
    
    for i = 1 : 2

      fc_rest = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),1,ifoi,:),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,1,ifoi,:),7),3);
      fc_task = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),2,ifoi,:),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,2,ifoi,:),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(2,2,i); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Blocks 1&2)'); else; title('Dpz - Pbo (Blocks 1&2)'); end
      
      colormap(cmap);
      tp_colorbar('Difference in correlation');
      
    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_f%d_v%d.pdf',ifoi,v))

  end
end

%% DRUG CONTRASTS FOR EACH BLOCK
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
cond = [2 3];

SUBJ = 1:28;
FOI = 6:7;
CLIM = [-0.02 0.02];

if v == 12
  
  [~,front_to_back] = sort(grid(:,2),'descend');
  
  left  = find(grid(:,1)<0);
  right = find(grid(:,1)>0);
  for ifoi = FOI
    
    figure; set(gcf,'color','w')
    
    for i = 1 : 2
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),1,ifoi,1),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,1,ifoi,1),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_rest = [fc1 fc2; fc3 fc4];
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),2,ifoi,1),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,2,ifoi,1),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_task = [fc1 fc2; fc3 fc4];
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      
      subplot(2,2,i); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #1)'); else; title('Dpz - Pbo (Block #1)'); end
      colormap(cmap)
      tp_colorbar('Correlation');
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),1,ifoi,2),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,1,ifoi,2),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_rest = [fc1 fc2; fc3 fc4];
      
      fc = nanmean(nanmean(cleandat(:,:,SUBJ,cond(i),2,ifoi,2),7),3)-nanmean(nanmean(cleandat(:,:,SUBJ,1,2,ifoi,2),7),3);
      
      fc1 = fc(front_to_back(left),front_to_back(left));
      fc2 = fc(front_to_back(left),front_to_back(right));
      fc3 = fc(front_to_back(right),front_to_back(left));
      fc4 = fc(front_to_back(right),front_to_back(right));
      
      fc_task = [fc1 fc2; fc3 fc4];
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(2,2,i+2); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #2)'); else; title('Dpz - Pbo (Block #2'); end
      
      colormap(cmap);
      tp_colorbar('Correlation');
      
    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_indivblocks_f%d_v%d.pdf',ifoi,v))

  end
else
  for ifoi = FOI
    figure; set(gcf,'color','w')
    
    for i = 1 : 2
      
      fc_rest = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),1,ifoi,1),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,1,ifoi,1),7),3);
      fc_task = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),2,ifoi,1),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,2,ifoi,1),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      
      subplot(2,2,i); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #1)'); else; title('Dpz - Pbo (Block #1)'); end
      colormap(cmap)
      tp_colorbar('Correlation');
      
      fc_rest = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),1,ifoi,2),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,1,ifoi,2),7),3);
      fc_task = nanmean(nanmean(cleandat(idx,idx,SUBJ,cond(i),2,ifoi,2),7),3)-nanmean(nanmean(cleandat(idx,idx,SUBJ,1,2,ifoi,2),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(2,2,i+2); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #2)'); else; title('Dpz - Pbo (Block #2'); end
      
      colormap(cmap);
      tp_colorbar('Correlation');
      
    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_indivblocks_f%d_v%d.pdf',ifoi,v))

  end
end



%% CORRELATIONS OF EFFECTS ACROSS BLOCKS AND CONDITIONS
ipharm = 3; icont = 1; clear r; clc

% pharma effect: correlation between blocks
for ifoi = 1 : 13
  for isubj = 1:28
    
    fc1 = cleandat(:,:,isubj,ipharm,icont,ifoi,1)-cleandat(:,:,isubj,1,icont,ifoi,1);
    fc2 = cleandat(:,:,isubj,ipharm,icont,ifoi,2)-cleandat(:,:,isubj,1,icont,ifoi,2);
    
    r(isubj,ifoi) = corr(fc1(mask),fc2(mask));
    
  end
end

fprintf('Mean correlation at 11 Hz: r = %.3f | p = %.3f\n',nanmean(r(:,6)))

% pharma effect: correlation between conditions
for ifoi = 1 : 13
  for isubj = 1:28
    
    fc1 = nanmean(cleandat(:,:,isubj,ipharm,1,ifoi,:),7)-nanmean(cleandat(:,:,isubj,1,1,ifoi,1),7);
    fc2 = nanmean(cleandat(:,:,isubj,ipharm,2,ifoi,:),7)-nanmean(cleandat(:,:,isubj,1,2,ifoi,1),7);
    
    r(isubj,ifoi) = corr(fc1(mask),fc2(mask));
    
  end
end

%% PLOT RAW FC ON CORTICAL SURAFCE

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

% close all

ifoi = 6; ipharm = 2; icont = 1; iblock = 1;

% cmap =  hot
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);
% cmap = cmap;

par_drug = nanmean(nanmean(nanmean(cleandat(:,:,:,ipharm,icont,ifoi,iblock),7),3));
par_plac = nanmean(nanmean(nanmean(cleandat(:,:,:,1,icont,ifoi,iblock),7),3));

[h,p]=ttest(nanmean(nanmean(cleandat(:,:,:,ipharm,icont,ifoi,iblock),7)),nanmean(nanmean(cleandat(:,:,:,1,icont,ifoi,iblock),7)),'dim',3);
par = 100*(par_drug-par_plac)./par_plac;
% par = par.*(p<0.4);

% % para = [];
% para.clim = [min(par) max(par)];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-15 15];
para.fn = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_cortexsurf_f%d_cont%d_v%d.pdf',ifoi,ipharm,v))

%% CORRELATIONS BETWEEN BLOCKS
mask = logical(tril(ones(400,400),-1));

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        dat_res(isubj,m,iblock) = 1;
      catch me
        %           warning(fprintf('S%d,M%d,B%d',isubj,m,iblock))
        dat_res(isubj,m,iblock) = 0;
      end
      try
        load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        dat_cnt(isubj,m,iblock) = 1;
      catch me
        %           warning(fprintf('S%d,M%d,B%d',isubj,m,iblock))
        dat_cnt(isubj,m,iblock) = 0;
      end
    end
  end
end

%% CORRELATIONS FC WITH BEHAVIOR
% clear nan_idx
%
% SUBJLIST1  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% para.str_behav = 'count';
% behav = pconn_read_behavioral_data(SUBJLIST1,para);
% behav_cnt = behav;
%
% para.str_behav = 'numb_switches';
% behav = pconn_read_behavioral_data(SUBJLIST1,para);
% behav_bttn = behav;
% behav_bttn = permute(behav_bttn,[2 1 3]);
%
% ifoi = 6; ipharm = 2; icont = 2; iblock = 1:2;
% cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);
%
% par_drug = squeeze(nanmean(nanmean(cleandat(:,:,:,ipharm,icont,ifoi,iblock),7)));
% par_plac = squeeze(nanmean(nanmean(cleandat(:,:,:,1,icont,ifoi,iblock),7)));
% par_diff = par_drug - par_plac;
% nan_idx1 = isnan(par_diff(1,:));
%
% par_behav = nanmean(behav_cnt(:,:,iblock),3)';
% par_behav = par_behav(:,ipharm)-par_behav(:,1);
% nan_idx2   = isnan(par_behav);
%
% nan_idx = ~(nan_idx1(:) |nan_idx2(:));
%
% par_behav = repmat(par_behav,[1 400])';
%
% [r,p] = corr(par_behav(:,nan_idx)',par_diff(:,nan_idx)'); r = r(1,:); p = p(1,:);
%
% para.cmap = cmap;
% para.grid = grid;
% para.dd = 0.75;
% para.clim = [-0.3 0.3];
% para.fn = sprintf('~/pupmod/plots/test.png');
% tp_plot_surface(r,para)
%


%% A FEW VERY BASIC ANALYSES
clear r

% Correlate ATX effect across blocks
eff1 = squeeze(fc(:,:,:,2,2,6,1)-fc(:,:,:,1,2,6,1));
eff2 = squeeze(fc(:,:,:,2,2,6,2)-fc(:,:,:,1,2,6,2));

for isubj = 1 : 28
  tmp1 = eff1(:,:,isubj);
  tmp2 = eff2(:,:,isubj);
  
  r(isubj) = corr(tmp1(mask),tmp2(mask));
  
end

[~,p]=ttest(r)
nanmean(r)

% averaged across vertices
eff1 = squeeze(nanmean(nanmean(squeeze(fc(:,:,:,2,2,6,1)-fc(:,:,:,1,2,6,1)),1),2));
eff2 = squeeze(nanmean(nanmean(squeeze(fc(:,:,:,2,2,6,2)-fc(:,:,:,1,2,6,2)),1),2));

eff1(27,:) = []; eff2(27,:) = [];

[r,p]=corr(eff1,eff2)

%% STH IS STRANGE WITH SUBJECT 13

figure; set(gcf,'color','w'); hold on
plot(squeeze(nanmean(nanmean(fc(:,:,13,1,1,:,1),1),2)))
plot(squeeze(nanmean(nanmean(fc(:,:,13,1,1,:,2),1),2)),'b--')

plot(squeeze(nanmean(nanmean(fc(:,:,13,1,2,:,1),1),2)),'r')
plot(squeeze(nanmean(nanmean(fc(:,:,13,1,2,:,2),1),2)),'r--')

plot(squeeze(nanmean(nanmean(fc(:,:,13,2,2,:,1),1),2)),'m')
plot(squeeze(nanmean(nanmean(fc(:,:,13,2,2,:,2),1),2)),'m--')

%%
clear r
mask      = logical(tril(ones(400,400),-1));
for ifoi = 1:13;

for i = 1 : 400
  i
  for j = 1 : 400
    
    r(i,j,ifoi) = corr(squeeze(fc(i,j,SUBJ,1,1,ifoi,1)),squeeze(fc(i,j,SUBJ,1,1,ifoi,2)));
    
  end
end
end
    %%
    
%     p = (sum((p_cnt2<0.05))+sum((p_cnt1<0.05)))./2;
para =[]
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.03 0.03];
para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_placebo_avg.png');
tp_plot_surface(nanmean(rr(:,:,8).*(pp(:,:,8)<0.05))',para)
%%
fc_mean   = squeeze(fc(:,:,:,:,2,ifoi,:));
fc_mean   = zscore(reshape(fc_mean,[400 400 28 6]),0,4);
fc_mean   = reshape(fc_mean,[400 400 28 3 2]);
fc_mean   = nanmean(fc_mean,5);
    
    
    
  

