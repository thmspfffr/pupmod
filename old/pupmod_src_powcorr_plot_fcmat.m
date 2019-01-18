%% pupmod_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

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
  
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v))


%%
  figure; set(gcf,'color','w')

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
  colormap(redblue)

end

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_powcorr_fcmat_f%d_v%d.pdf',ifoi,v))
%%
