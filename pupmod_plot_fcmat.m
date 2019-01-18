%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

% Last changed: 13-11-2018

clear

v = 1;

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath(genpath('/home/gnolte/meth'));

load sa_meg_template;

grid  = select_chans(sa_meg_template.grid_cortex3000,400);
  
cleandat = pupmod_loadpowcorr(v,1);

%%

  figure; set(gcf,'color','w')

  load redblue.mat
  
cmap = cbrewer('div', 'RdBu', 256,'pchip');
cmap = cmap(end:-1:1,:);  
  
ifoi = 7;
cond = [2 3]; 
[~,front_to_back] = sort(grid(:,2),'descend');

left  = find(grid(:,1)<0);
right = find(grid(:,1)>0);

for i = 1 : 2

  fc = 100*(nanmean(cleandat(:,:,:,cond(i),1,ifoi),3)-nanmean(cleandat(:,:,:,1,1,ifoi),3))./nanmean(cleandat(:,:,:,1,1,ifoi),3);
  
  fc1 = fc(front_to_back(left),front_to_back(left));
  % fc1 = zeros(201,201);

  fc2 = fc(front_to_back(left),front_to_back(right));
  % fc2 = 10*ones(201,199);

  fc3 = fc(front_to_back(right),front_to_back(left));
  fc4 = fc(front_to_back(right),front_to_back(right));

  fc_rest = [fc1 fc2; fc3 fc4];

  fc = 100*(nanmean(cleandat(:,:,:,cond(i),2,ifoi),3)-nanmean(cleandat(:,:,:,1,2,ifoi),3))./nanmean(cleandat(:,:,:,1,2,ifoi),3);

  fc1 = fc(front_to_back(left),front_to_back(left));
  % fc1 = zeros(201,201);
  fc2 = fc(front_to_back(left),front_to_back(right));
  % fc2 = 10*ones(201,199);
  fc3 = fc(front_to_back(right),front_to_back(left));
  fc4 = fc(front_to_back(right),front_to_back(right));

  fc_task = [fc1 fc2; fc3 fc4];

  fc = tril(fc_rest,-1)+triu(fc_task,1);
  % subplot(1,2,1); imagesc(fc_rest,[-30 30]); axis square off
  subplot(1,2,i); imagesc(fc,[-30 30]); axis square off
  colormap(cmap)

end

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_f%d_v%d.pdf',ifoi,v))
%%

addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load redblue.mat
  
cmap = cbrewer('div', 'RdBu', 256,'pchip');
cmap = cmap(end:-1:1,:);  
  
ifoi = 6;
cond = [2 3]; 

fc = nanmean(cleandat(:,:,:,cond(i),1,ifoi),3);

para = [];
para.n = 90;
para.transfer = 'to_bcn';

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

[fc, ~,lab] = tp_match_aal(para,fc); 
fc = fc(include_bcn,include_bcn);
lab = lab(include_bcn);

figure; set(gcf,'color','w')
set(gcf,'Position',[100 100 1000 1000])
  
cmap = cbrewer('div', 'RdBu', 256,'pchip');
cmap = cmap(end:-1:1,:);  

imagesc(fc,[0 0.10]); colormap(plasma)
set(gca,'xtick',1:76,'xticklabel',erase(lab,'_'),'tickdir','out')
xtickangle(90)
set(gca,'ytick',1:76,'yticklabel',erase(lab,'_'),'tickdir','out')
axis square 
tp_editplots

print(gcf,'-dpdf',sprintf('~/pupmod/plots/clopupmod_powcorr_fcmat_aal_f%d_v%d.pdf',ifoi,v))

