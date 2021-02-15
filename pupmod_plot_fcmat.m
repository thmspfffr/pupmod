%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

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

fc = pupmod_loadpowcorr(v,SUBJLIST,1);

%% NEW PLOTS (30-04-2020)


% PLOT FC MATRIX FOR DIFFERENT CONDITIONS
% significant v23: atx 10/11/12 & 18, dpz 13/14
ifoi = 9;

set(0,'defaultfigurerenderer','painters')
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

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
  
  fc_rest_tmp = [fc1 fc2; fc3 fc4];
  fc_rest_tmp = triu(fc_rest_tmp,1);
  fc_rest(:,:,i) = rot90(flipud(fc_rest_tmp),-1);
  
  figure_w;
  imagesc(fc_rest(:,:,i),[0.03 0.10]); axis square off
  colormap(plasma)
  
  print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_rest_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))
  
  
  fc_tmp = nanmean(nanmean(nanmean(fc(:,:,:,cond(i),2,ifoi,:),6),7),3);
  %
  fc1 = fc_tmp(front_to_back(left),front_to_back(left));
  fc2 = fc_tmp(front_to_back(left),front_to_back(right));
  fc3 = fc_tmp(front_to_back(right),front_to_back(left));
  fc4 = fc_tmp(front_to_back(right),front_to_back(right));
  
  fc_task_tmp = [fc1 fc2; fc3 fc4];
  fc_task_tmp = triu(fc_task_tmp,1);
  fc_task(:,:,i) = fc_task_tmp;
  
  figure_w;
  imagesc(fc_task(:,:,i) ,[0.03 0.10]); axis square off
  colormap(plasma)
  
  print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_task_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))
  
  
end

cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);

figure_w;
imagesc(fc_rest(:,:,2)-fc_rest(:,:,1),[-0.02 0.02]); axis square off
colormap(cmap)

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_diff_atx_rest_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))

figure_w;
imagesc(fc_rest(:,:,3)-fc_rest(:,:,1),[-0.02 0.02]); axis square off
colormap(cmap)

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_diff_dpz_rest_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))

figure_w;
imagesc(fc_task(:,:,2)-fc_task(:,:,1),[-0.02 0.02]); axis square off

colormap(cmap)

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_diff_atx_task_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))

figure_w;
imagesc(fc_task(:,:,3)-fc_task(:,:,1),[-0.02 0.02]); axis square off
colormap(cmap)

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_raw_fcmat_diff_dpz_task_f%s_c%d_v%d.eps',regexprep(num2str(ifoi),' ',''),cond(i),v))

    

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

%% CORRELATE DRUG EFFECTS REST AND TASK

for isubj = 1 : 28
  for ifoi = 1 : 25
    
    atx_rest = squeeze(nanmean(fc(:,:,isubj,2,1,ifoi),2)-nanmean(fc(:,:,isubj,1,1,ifoi),2));
    atx_task = squeeze(nanmean(fc(:,:,isubj,2,2,ifoi),2)-nanmean(fc(:,:,isubj,1,2,ifoi),2));
    
    dpz_rest = squeeze(nanmean(fc(:,:,isubj,3,1,ifoi),2)-nanmean(fc(:,:,isubj,1,1,ifoi),2));
    dpz_task = squeeze(nanmean(fc(:,:,isubj,3,2,ifoi),2)-nanmean(fc(:,:,isubj,1,2,ifoi),2));
    
    [corr_taskrest_atx(isubj,ifoi), p_atx(isubj,ifoi)] = corr(atx_rest,atx_task);
    [corr_taskrest_dpz(isubj,ifoi), p_dpz(isubj,ifoi)]= corr(dpz_rest,dpz_task);
  
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

ifoi = 5:9; ipharm = 2; icont = 2; iblock = 1;

% cmap =  hot
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);
% cmap = cmap;

par_drug = nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,icont,ifoi,iblock),7),3),6));
par_plac = nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,icont,ifoi,iblock),7),3),6));

[h,p]=ttest(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,icont,ifoi,iblock),7),6)),nanmean(nanmean(nanmean(fc(:,:,:,1,icont,ifoi,iblock),7),6)),'dim',3);
par = 100*(par_drug-par_plac)./par_plac;
  par = par.*(p<0.05);

% % para = [];
% para.clim = [min(par) max(par)];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.clim = [-20 20];
para.fn = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)

% print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_cortexsurf_f%d_cont%d_v%d.pdf',ifoi,ipharm,v))

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

%% PLOT FOR METHODS FIGURE


figure; set(gcf,'color','w')
cond = [1 2 3];
[~,front_to_back] = sort(grid(:,2),'descend');
left  = find(grid(:,1)<0);
right = find(grid(:,1)>0);

isubj = 18;
icond = 1;

  
fc_tmp = nanmean(nanmean(nanmean(fc(:,:,isubj,1,icond,ifoi,:),6),7),3);

fc1 = fc_tmp(front_to_back(left),front_to_back(left));
fc2 = fc_tmp(front_to_back(left),front_to_back(right));
fc3 = fc_tmp(front_to_back(right),front_to_back(left));
fc4 = fc_tmp(front_to_back(right),front_to_back(right));

fc_rest_tmp = [fc1 fc2; fc3 fc4];
% fc_rest(:,:,i) = triu(fc_rest_tmp,1);

subplot(1,3,i); imagesc(fc_rest_tmp,[0.03 0.10]); axis square off
colormap(plasma)
  
 

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_fcmat_methodsfig_isubj%d_cond%d_v%d.eps',isubj,icond,v))


figure; set(gcf,'color','w')
cond = [1 2 3];
[~,front_to_back] = sort(grid(:,2),'descend');
left  = find(grid(:,1)<0);
right = find(grid(:,1)>0);

for isubj = 1 : 28
  icond = 1;
ifoi = 13
  
fc_tmp = nanmean(nanmean(nanmean(fc(:,:,isubj,1,2,ifoi,:),6),7),3)-nanmean(nanmean(nanmean(fc(:,:,isubj,1,1,ifoi,:),6),7),3);

fc1 = fc_tmp(front_to_back(left),front_to_back(left));
fc2 = fc_tmp(front_to_back(left),front_to_back(right));
fc3 = fc_tmp(front_to_back(right),front_to_back(left));
fc4 = fc_tmp(front_to_back(right),front_to_back(right));

fc_rest_tmp = [fc1 fc2; fc3 fc4];

subplot(1,3,i); imagesc(fc_rest_tmp,[-0.1 0.1]); axis square off
colormap(cmap)
  
print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_fcmat_methodsfig_isubj%d_cond%d_v%d.eps',isubj,icond,v))

end
%%
    
