%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

% Last changed: 13-11-2018

clear

v = 25;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

fprintf('Loading grid...\n')
grid  = select_chans(sa_meg_template.grid_cortex3000,400); fprintf('Loading grid... Done\n')
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

fc = pupmod_loadpowcorr(v,SUBJLIST,0);

if v == 20
  tmp = tp_create_grid('vtpm');
  clim = [-0.025 0.025];
  idx = 1:46;
  iismember(idx,[21 22 23 44 45 46]);
  reg = tmp.tissuelabel_4mm(idx);
end
%% PLOT FC MATRIX FOR DIFFERENT CONDITIONS
% significant v23: atx 10/11/12 & 18, dpz 13/14
ifoi = 13:14;


cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

if v == 12 | v == 19 | v==23
%   for ifoi = FOI

    figure; set(gcf,'color','w')
    cond = [1 2 3];
    [~,front_to_back] = sort(grid(:,2),'descend');
    left  = find(grid(:,1)<0);
    right = find(grid(:,1)>0);
    for i = 1 : 3
      
      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,cond(i),1,ifoi,:),6),7),3);

      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));

      fc_rest = [fc1 fc2; fc3 fc4];

      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,cond(i),2,ifoi,:),6),7),3);

      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));

      fc_task = [fc1 fc2; fc3 fc4];

      fc_tmp = tril(fc_rest,-1)+triu(fc_task,1);
      
      min_val = min([min(fc_tmp(logical(triu(ones(400,400),1)))) min(fc_tmp(logical(tril(ones(400,400),-1))))]);
      max_val = max([max(fc_tmp(logical(triu(ones(400,400),1)))) max(fc_tmp(logical(tril(ones(400,400),-1))))]);
      
      
      subplot(1,3,i); imagesc(fc_tmp,[2*min_val 0.75*max_val]); axis square off
      colormap(plasma)

    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_f%s_v%d.pdf',regexprep(num2str(ifoi),' ',''),v))

    figure; set(gcf,'color','w');
    
    % PLOT DRUG CONTRASTS 
    % ----------------
    % Atomoxetine first
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,2,1,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,1,1,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_rest = [fc1 fc2; fc3 fc4];
    
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,2,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,1,2,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_task = [fc1 fc2; fc3 fc4];
    
    fc_tmp = tril(fc_rest,-1)+triu(fc_task,1);
    
%     min_val = min([min(fc_tmp(logical(triu(ones(400,400),1)))) min(fc_tmp(logical(tril(ones(400,400),-1))))]);
%     max_val = max([max(fc_tmp(logical(triu(ones(400,400),1)))) max(fc_tmp(logical(tril(ones(400,400),-1))))]);
    min_val = -0.02; max_val = 0.02;

    
    subplot(1,3,2); imagesc(fc_tmp,[min_val max_val]); axis square off
    colormap(cmap)
    
    % Donepezil second
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,3,1,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,1,1,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_rest = [fc1 fc2; fc3 fc4];
    
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,3,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,1,2,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_task = [fc1 fc2; fc3 fc4];
    
    fc_tmp = tril(fc_rest,-1)+triu(fc_task,1);
    
%     min_val = min([min(fc_tmp(logical(triu(ones(400,400),1)))) min(fc_tmp(logical(tril(ones(400,400),-1))))]);
%     max_val = max([max(fc_tmp(logical(triu(ones(400,400),1)))) max(fc_tmp(logical(tril(ones(400,400),-1))))]);
%     min_val = -0.03; max_val = 0.03;


    subplot(1,3,3); imagesc(fc_tmp,[min_val max_val]); axis square off
    colormap(cmap)
    
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_diff_fcmat_f%s_v%d.pdf',regexprep(num2str(ifoi),' ',''),v))

    % -------------
    % TASK VS REST
    % --------------------------------------------------
     % PLOT TASK VS REST CONTRASTS 
    % ----------------
    % Placebo first
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,1,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,1,1,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_rest = [fc1 fc2; fc3 fc4];
    
    subplot(1,3,1); imagesc(fc_rest,[min_val max_val]);axis square off
    colormap(cmap)
    
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,2,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,2,1,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_rest = [fc1 fc2; fc3 fc4];
    
    
    subplot(1,3,2); imagesc(fc_rest,[min_val max_val]);axis square off
    colormap(cmap)
    
    fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,3,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,:,3,1,ifoi,:),6),7),3);

    fc1 = fc_tmp(front_to_back(left),front_to_back(left));
    fc2 = fc_tmp(front_to_back(left),front_to_back(right));
    fc3 = fc_tmp(front_to_back(right),front_to_back(left));
    fc4 = fc_tmp(front_to_back(right),front_to_back(right));

    fc_rest = [fc1 fc2; fc3 fc4];
    
    
    subplot(1,3,3); imagesc(fc_rest,[min_val max_val]);axis square off
    colormap(cmap)
    
elseif v== 20
  for ifoi = FOI

    figure; set(gcf,'color','w')
    cond = [1 2 3];
    
    for i = 1 : 3
      
      fc_rest = nanmean(nanmean(fc(idx,idx,:,cond(i),1,ifoi,:),7),3);
      fc_task = nanmean(nanmean(fc(idx,idx,:,cond(i),2,ifoi,:),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
%       if ifoi ~= 7
      subplot(1,3,3); imagesc(fc,[0.04 0.12]); axis square %off
      
      set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
      set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
      colormap(inferno); tp_editplots; set(gca,'FontSize',5)

    end
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_f%s_v%d.pdf',regexprep(num2str(ifoi),' ',''),v))
  end
end

%% DRUG CONTRASTS FOR EACH BLOCK
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
cond = [2 3];

SUBJ = 1:28;
ifoi = 13:14;
CLIM = [-0.02 0.02];

if v == 12 | v == 19 | v == 23
  
  [~,front_to_back] = sort(grid(:,2),'descend');
  
  left  = find(grid(:,1)<0);
  right = find(grid(:,1)>0);
%   for ifoi = FOI
    
    figure; set(gcf,'color','w')
    
    for i = 1 : 2
      
      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,SUBJ,cond(i),1,ifoi,1),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,1,ifoi,1),6),7),3);
      
      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));
      
      fc_rest = [fc1 fc2; fc3 fc4];
      
      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,SUBJ,cond(i),2,ifoi,1),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,2,ifoi,1),6),7),3);
      
      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));
      
      fc_task = [fc1 fc2; fc3 fc4];
      
      fc_tmp = tril(fc_rest,-1)+triu(fc_task,1);
      
      subplot(2,2,i); imagesc(fc_tmp,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #1)'); else; title('Dpz - Pbo (Block #1)'); end
      colormap(cmap)
      tp_colorbar('Correlation');
      
      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,SUBJ,cond(i),1,ifoi,2),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,1,ifoi,2),6),7),3);
      
      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));
      
      fc_rest = [fc1 fc2; fc3 fc4];
      
      fc_tmp = nanmean(nanmean(nanmean(fc(:,:,SUBJ,cond(i),2,ifoi,2),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,SUBJ,1,2,ifoi,2),6),7),3);
      
      fc1 = fc_tmp(front_to_back(left),front_to_back(left));
      fc2 = fc_tmp(front_to_back(left),front_to_back(right));
      fc3 = fc_tmp(front_to_back(right),front_to_back(left));
      fc4 = fc_tmp(front_to_back(right),front_to_back(right));
      
      fc_task = [fc1 fc2; fc3 fc4];
      
      fc_tmp = tril(fc_rest,-1)+triu(fc_task,1);
      subplot(2,2,i+2); imagesc(fc_tmp,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #2)'); else; title('Dpz - Pbo (Block #2'); end
      
      colormap(cmap);
      tp_colorbar('Correlation');
      
%     end

    end
      print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_indivblocks_f%s_v%d.pdf',regexprep(num2str(ifoi),' ',''),v))

else
  for ifoi = FOI
    figure; set(gcf,'color','w')
    
    for i = 1 : 2
      
      fc_rest = nanmean(nanmean(fc(idx,idx,SUBJ,cond(i),1,ifoi,1),7),3)-nanmean(nanmean(fc(idx,idx,SUBJ,1,1,ifoi,1),7),3);
      fc_task = nanmean(nanmean(fc(idx,idx,SUBJ,cond(i),2,ifoi,1),7),3)-nanmean(nanmean(fc(idx,idx,SUBJ,1,2,ifoi,1),7),3);
      
      fc = tril(fc_rest,-1)+triu(fc_task,1);
      
      subplot(2,2,i); imagesc(fc,CLIM); axis square off;
      if i == 1; title('Atx - Pbo (Block #1)'); else; title('Dpz - Pbo (Block #1)'); end
      colormap(cmap)
      tp_colorbar('Correlation');
      
      fc_rest = nanmean(nanmean(fc(idx,idx,SUBJ,cond(i),1,ifoi,2),7),3)-nanmean(nanmean(fc(idx,idx,SUBJ,1,1,ifoi,2),7),3);
      fc_task = nanmean(nanmean(fc(idx,idx,SUBJ,cond(i),2,ifoi,2),7),3)-nanmean(nanmean(fc(idx,idx,SUBJ,1,2,ifoi,2),7),3);
      
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
    
    fc1 = fc(:,:,isubj,ipharm,icont,ifoi,1)-fc(:,:,isubj,1,icont,ifoi,1);
    fc2 = fc(:,:,isubj,ipharm,icont,ifoi,2)-fc(:,:,isubj,1,icont,ifoi,2);
    
    r(isubj,ifoi) = corr(fc1(mask),fc2(mask));
    
  end
end

fprintf('Mean correlation at 11 Hz: r = %.3f | p = %.3f\n',nanmean(r(:,6)))

% pharma effect: correlation between conditions
for ifoi = 1 : 13
  for isubj = 1:28
    
    fc1 = nanmean(fc(:,:,isubj,ipharm,1,ifoi,:),7)-nanmean(fc(:,:,isubj,1,1,ifoi,1),7);
    fc2 = nanmean(fc(:,:,isubj,ipharm,2,ifoi,:),7)-nanmean(fc(:,:,isubj,1,2,ifoi,1),7);
    
    r(isubj,ifoi) = corr(fc1(mask),fc2(mask));
    
  end
end

%% PLOT RAW FC ON CORTICAL SURAFCE

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v9.mat
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

par_drug = nanmean(nanmean(nanmean(fc(:,:,:,ipharm,icont,ifoi,iblock),7),3));
par_plac = nanmean(nanmean(nanmean(fc(:,:,:,1,icont,ifoi,iblock),7),3));

[h,p]=ttest(nanmean(nanmean(fc(:,:,:,ipharm,icont,ifoi,iblock),7)),nanmean(nanmean(fc(:,:,:,1,icont,ifoi,iblock),7)),'dim',3);
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
        load(sprintf('~/pp/proc/pp_sens_fc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        dat_res(isubj,m,iblock) = 1;
      catch me
        %           warning(fprintf('S%d,M%d,B%d',isubj,m,iblock))
        dat_res(isubj,m,iblock) = 0;
      end
      try
        load(sprintf('~/pp/proc/pp_cnt_sens_fc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        dat_cnt(isubj,m,iblock) = 1;
      catch me
        %           warning(fprintf('S%d,M%d,B%d',isubj,m,iblock))
        dat_cnt(isubj,m,iblock) = 0;
      end
    end
  end
end


%% CORRELATIONS ACROSS FREQS
mask    = logical(triu(ones(400,400),1));
for icond = 1 : 1
for ifreq1 = 1 : 25
  ifreq1
  for ifreq2 = 1 : 25
    
    tmp1 = nanmean(cleandat(:,:,:,2,icond,ifreq1),3); tmp1 = tmp1(mask);
    tmp2 = nanmean(cleandat(:,:,:,2,icond,ifreq2),3); tmp2 = tmp2(mask);
    
    r(ifreq1,ifreq2,icond) = corr(tmp1,tmp2);
        
  end
end
end

%%


% A FEW VERY BASIC ANALYSES
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
for ifoi = 1:13
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
  
%% TASK VS REST
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

ipharm = 3;
SUBJ = 1:28;

% for ipharm = 1 : 3
% 
%   for ifoi = 1 : 13
% ifoi
% d_fc1       = (squeeze(nanmean(fc(:,:,SUBJ,1,2,ifoi,:),7))-squeeze(nanmean(fc(:,:,SUBJ,1,1,ifoi,:),7)))./squeeze(nanmean(fc(:,:,SUBJ,1,1,ifoi,:),7));%./squeeze(nanmean(fc(:,:,SUBJ,1,1,ifoi,:),7));
% d_fc2       = (squeeze(nanmean(fc(:,:,SUBJ,2,2,ifoi,:),7))-squeeze(nanmean(fc(:,:,SUBJ,2,1,ifoi,:),7)));%./squeeze(nanmean(fc(:,:,SUBJ,2,1,ifoi,:),7));
% d_fc3       = (squeeze(nanmean(fc(:,:,SUBJ,3,2,ifoi,:),7))-squeeze(nanmean(fc(:,:,SUBJ,3,1,ifoi,:),7)));%./squeeze(nanmean(fc(:,:,SUBJ,3,1,ifoi,:),7));

d_fc1=(nanmean(fc(:,:,SUBJ,1,2,ifoi,:),7)-nanmean(fc(:,:,SUBJ,1,1,ifoi,:),7));

% [h,~,~,s]=ttest(squeeze(d_fc2),squeeze(d_fc1),'dim',3,'alpha',0.01)

%   imagesc(nanmean(d_fc1,3),[-3 3])
%   colormap(cmap)
% 
%   [h,~,~,s]=ttest(d_fc,zeros(400,400,28),'dim',3);
% 
%   pos(ifoi,ipharm) = sum(h(mask)>0 & s.tstat(mask)>0)/sum(mask(:));
%   neg(ifoi,ipharm) = sum(h(mask)>0 & s.tstat(mask)<0)/sum(mask(:));
% 
% 
%   end
% end

%%
clear pow18_res
% isubj = 6; m = 3; 
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    load(sprintf('/home/tpfeffer/pp/proc/conn/pp_task_src_powcorr_test_s%d_m%d_v19.mat',isubj,im))
    pow18(:,:,:,isubj,m,:) = single(powcorr); clear powcorr
    load(sprintf('/home/tpfeffer/pp/proc/conn/pp_src_powcorr_test_s%d_m%d_v19.mat',isubj,im))
    pow18_res(:,:,:,isubj,m,:) = single(powcorr); clear powcorr
%     load(sprintf('/home/tpfeffer/pupmod/proc/conn/pupmod_task_src_powcorr_s%d_m%d_b1_f6_v12.mat',isubj,im))
%     pow1(:,:,1) = powcorr; clear powcorr
%     load(sprintf('/home/tpfeffer/pupmod/proc/conn/pupmod_task_src_powcorr_s%d_m%d_b2_f6_v12.mat',isubj,im))
%     p1(:,:,2) = powcorr; clear powcorr
%     pow1(:,:,isubj,m,:) = p1;
%   for iblock = 1 : 2
%     for ifoi = 1: 13
%     load(sprintf('/home/tpfeffer/pupmod/proc/conn/pupmod_src_powcorr_s%d_m%d_b%d_f%d_v12.mat',isubj,im,iblock,ifoi))
%     pow12_res(:,:,iblock,isubj,m,ifoi) = single(powcorr); clear powcorr
%     end
%   end
  end
end

pow18 = pow18(:,:,:,SUBJLIST,:,:);
pow18_res = pow18_res(:,:,:,SUBJLIST,:,:);
% pow12_res = pow12_res(:,:,:,SUBJLIST,:,:);

% powres = powres(:,:,:,SUBJLIST,:,:);
%% ALTERED CORR

iblock = 2;

for ifoi = 1 :25
  ifoi
  [h,~,~,s]=ttest(squeeze(fc(:,:,:,2,2,ifoi,iblock)),squeeze(fc(:,:,:,1,2,ifoi,iblock)),'dim',3);
  pos_atx(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_atx(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  [h,~,~,s]=ttest(squeeze(fc(:,:,:,3,2,ifoi,iblock)),squeeze(fc(:,:,:,1,2,ifoi,iblock)),'dim',3);
  pos_dpz(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_dpz(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
%   
  [h,~,~,s]=ttest(squeeze(fc(:,:,:,2,1,ifoi,iblock)),squeeze(fc(:,:,:,1,1,ifoi,iblock)),'dim',3);
  pos_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  [h,~,~,s]=ttest(squeeze(fc(:,:,:,3,1,ifoi,iblock)),squeeze(fc(:,:,:,1,1,ifoi,iblock)),'dim',3);
  pos_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  
%   [h,~,~,s]=ttest(squeeze(nanmean(pow12_res(:,:,:,:,2,ifoi),3)),squeeze(nanmean(pow12_res(:,:,:,:,1,ifoi),3)),'dim',3);
%   pos_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
%   neg_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
%   [h,~,~,s]=ttest(squeeze(nanmean(pow12_res(:,:,:,:,3,ifoi),3)),squeeze(nanmean(pow12_res(:,:,:,:,1,ifoi),3)),'dim',3);
%   pos_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
%   neg_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
end


%%
isubj = 1:28; iblock = 1:2; m = 2; foi =16;

cond = 2;

tmp_fc(:,:,:,:,1) = squeeze(nanmean(nanmean(pow18_res(:,:,iblock,isubj,:,foi),3),6));
tmp_fc(:,:,:,:,2) = squeeze(nanmean(nanmean(pow18(:,:,iblock,isubj,:,foi),3),6));

d_atx = squeeze(nanmean(tmp_fc(:,:,:,2,:)-tmp_fc(:,:,:,1,:),3));
d_dpz = squeeze(nanmean(tmp_fc(:,:,:,3,:)-tmp_fc(:,:,:,1,:),3));


para =[];
% para.cmap = plasma;
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.01 0.01]
% para.clim = [min(nanmean(d18,2)) max(nanmean(d18,2))]
para.fn = '~/test.png';
% para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
% tp_plot_surface(nanmean(d_dpz(:,:,cond),2),para);

% para.clim = [min(nanmean(d17,2)) max(nanmean(d17,2))]
% para.fn = '~/test.png';
% % para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
% tp_plot_surface(nanmean(d17,2),para);

% para.clim = [min(nanmean(d3,2)) max(nanmean(d3,2))]
para.fn = '~/test.png';
% para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
% tp_plot_surface(nanmean(d3,2),para);

cond = 2; ipharm = 1;
d_fc = squeeze(nanmean(tmp_fc(:,:,:,ipharm,cond),2));
d_behav = repmat(nanmean(behav_cnt(ipharm,:,:),3),[400 1]);

[r,p]=corr(d_behav',d_fc'); r=r(1,:);p=p(1,:);

para =[];
% para.cmap = plasma;
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.2 0.2]
% para.clim = [min(nanmean(d18,2)) max(nanmean(d18,2))]
para.fn = '~/test.png';
tp_plot_surface(r,para);


%%
% m = 3; foi = ;
% fc = 
foi =1
for isubj = 1:28
  isubj
%   for m = 1:3
%     for foi = 1 : 9
  

%   t1 = fc(:,:,find(ismember(SUBJLIST1,SUBJLIST(isubj))),2,2,6,1)-fc(:,:,find(ismember(SUBJLIST1,SUBJLIST(isubj))),1,2,6,1);t1 = t1(mask);
%   t2 = fc(:,:,find(ismember(SUBJLIST1,SUBJLIST(isubj))),2,2,6,2)-fc(:,:,find(ismember(SUBJLIST1,SUBJLIST(isubj))),1,2,6,2); t2 = t2(mask);
% % t1 = t1(mask);
% t2 = t2(mask);
  %   m1(isubj) = squeeze(nanmean(nanmean(pow1(:,:,isubj,m,1),2),1));
%   m2(isubj) = squeeze(nanmean(nanmean(pow1(:,:,isubj,m,2),2),1));
  
%   r1(isubj) = corr(t1,t2);
  
  t1 = nanmean(pow2(:,:,1,isubj,2,foi),6);%-nanmean(pow2(:,:,1,isubj,1,foi),6); t1 = t1(mask);
  t2 = nanmean(pow2(:,:,2,isubj,2,foi),6);%-nanmean(pow2(:,:,2,isubj,1,foi),6); t2 = t2(mask);
  t1 = t1(mask);
t2 = t2(mask);
  r2(isubj) = corr(t1,t2);
  
%   t1 = nanmean(pow13(:,:,1,isubj,2,1:9),6);%-nanmean(pow13(:,:,1,isubj,1,foi),6); t1 = t1(mask);
%   t2 = nanmean(pow13(:,:,2,isubj,2,1:9),6);%-nanmean(pow13(:,:,2,isubj,1,foi),6); t2 = t2(mask);
%   t1 = t1(mask);
% t2 = t2(mask);
%   r3(isubj) = corr(t1,t2);
%   end
% end
end
  
  

%%
 isubj = 1:27;
 foi = 9:11;
 iblock = 1:2;
% d1 = nanmean(nanmean(fc(:,:,find(ismember(SUBJLIST1,SUBJLIST)),2,2,6,iblock),7)-nanmean(fc(:,:,find(ismember(SUBJLIST1,SUBJLIST)),1,2,6,iblock),7),3);
d2 = nanmean(nanmean(nanmean(nanmean(pow18(:,:,iblock,isubj,2,foi),3),4),6)-nanmean(nanmean(nanmean(pow18(:,:,iblock,isubj,1,foi),3),4),6),3);
% d3 = nanmean(nanmean(nanmean(nanmean(powres(:,:,iblock,isubj,2,foi),3),4),6)-nanmean(nanmean(nanmean(powres(:,:,iblock,isubj,1,foi),3),4),6),3);

% 
% [h1,p1] = ttest(nanmean(nanmean(pow1(:,:,isubj,2,:),1),5),nanmean(nanmean(pow1(:,:,isubj,1,:),1),5),'dim',3);
% [h2,p2] = ttest(nanmean(nanmean(pow17(:,:,:,isubj,2,foi),3),6),nanmean(nanmean(pow2(:,:,:,isubj,1,foi),3),6),'dim',4);
% [h2,p2] = ttest(nanmean(nanmean(nanmean(pow2(:,:,:,isubj,2,foi),1),3),6),nanmean(nanmean(nanmean(pow2(:,:,:,isubj,1,foi),1),3),6),'dim',4);
% [hres,pres] = ttest(nanmean(nanmean(powres(:,:,:,isubj,2,foi),3),6),nanmean(nanmean(powres(:,:,:,isubj,1,foi),3),6),'dim',4);

% figure;
% subplot(1,2,1)
% imagesc(d1,[-0.02 0.02]);
% subplot(1,2,2)
% imagesc(d2,[-0.02 0.02]); colormap(cmap)
%
% corr(nanmean(d1,2),nanmean(d2,2))
para =[];
% para.cmap = plasma;
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.025 0.025];
% para.clim = [min(nanmean(d1,2)) max(nanmean(d1,2))]
para.fn = '~/test.png';
% para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(nanmean(d2,2),para);
% 
% 
% para.clim = [min(nanmean(d2,2)) max(nanmean(d2,2))]
para.fn = '~/test.png';
% para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
% tp_plot_surface(nanmean(d3,2),para);
% 
%%

para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt = behav;

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);
%%
foi = 9:11;
m = 3;

d_meg = squeeze(nanmean(nanmean(pow18(:,:,:,1:27,m,foi),6),3));%-squeeze(nanmean(nanmean(pow2(:,:,:,:,1,foi),6),3));
d_behav = nanmean(behav_cnt(m,find(ismember(SUBJLIST1,SUBJLIST)),:),3);%- nanmean(behav_cnt(1,find(ismember(SUBJLIST1,SUBJLIST)),:),3);
d_behav = permute(repmat(d_behav(:),[1 400 400]),[2 3 1]);

[r,p]=tp_corr(d_meg,d_behav,3);

para =[];
% para.cmap = plasma;
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-0.2 0.2]
% para.clim = [min(nanmean(d1,2)) max(nanmean(d1,2))]
para.fn = '~/test.png';
% para.fn = sprintf('~/pupmod/plots/pupmod_behav_fc_corr_with_behav_surf_f%d_b%s_pharm%s_v%d.png',ifoi,regexprep(num2str(iblock),' ',''),regexprep(num2str(ipharm),' ',''),v);
tp_plot_surface(nanmean(r),para);

%%
for isubj = 1:28
d_meg = squeeze(nanmean(nanmean(pow2(:,:,:,isubj,2,foi),6),3))-squeeze(nanmean(nanmean(pow2(:,:,:,isubj,1,foi),6),3));
d_meg = d_meg(mask);
d(isubj) = mean(d_meg(p2(mask)<1));

end


%% MIN MAX

% d1 = min(reshape(nanmean(pow1(:,:,:,2,1),5),[400 * 400 28]))-min(reshape(nanmean(pow1(:,:,:,2,2),5),[400 * 400 28])); 
% d2 = min(reshape(nanmean(pow2(:,:,:,2,1),5),[400 * 400 28]))-min(reshape(nanmean(pow2(:,:,:,2,2),5),[400 * 400 28])); 

d1 = max(reshape(nanmean(pow1(:,:,:,2,1),5),[400 * 400 28]))-max(reshape(nanmean(pow1(:,:,:,2,2),5),[400 * 400 28])); 
d2 = max(reshape(nanmean(pow2(:,:,:,2,1),5),[400 * 400 28]))-max(reshape(nanmean(pow2(:,:,:,2,2),5),[400 * 400 28])); 

range1 = squeeze(max(reshape(pow1(:,:,:,2,:),[400 * 400 28 2]))-min(reshape(pow1(:,:,:,2,:),[400 * 400 28 2])));
range2 = squeeze(max(reshape(pow2(:,:,:,2,:),[400 * 400 28 2]))-min(reshape(pow1(:,:,:,2,:),[400 * 400 28 2])));
