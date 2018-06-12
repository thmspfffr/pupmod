%% pupmod_all_src_genemaps
% COMPARES CHANGES IN FC TO GENE EXPRESSION MAPS PROVIDED BY RUDY

clear
v         = 12; % 400 voxels LCMV
  

outdir = '~/pupmod/proc/conn/';

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
load ~/pmod/matlab/gene_values.mat
load ~/pconn/proc/src/pconn_sa_s20_m1_b1_v10.mat
grid = sa.grid_cortex_lowres; clear sa
orig = [46 64 37];

% Transform gene maps to MNI
locs=2*(locs-repmat(orig,[size(locs,1) 1]))/10;

for igrid = 1 : size(grid,1)

  pos = grid(igrid,:);
  
  for ilocs = 1 : size(locs,1)
    d(ilocs) = sqrt( (locs(ilocs,1)-pos(1))^2 + (locs(ilocs,2)-pos(2))^2 + (locs(ilocs,3)-pos(3))^2 );
  end
  
  dist = d(d < 2);
  idx  = find(d<2);
  if isempty(idx) 
    w_gdat(igrid,:) = nan([1 size(gdat,2)]);
    continue
  end
  
  w = exp((-(dist).^2)./4)%/sum(exp((-(dist).^2)./4));
  
  w_gdat(igrid,:) = sum(gdat(idx,:) .* repmat(w(:),[1 size(gdat,2)]));
  

end

idx = ~isnan(w_gdat(:,1));
fc = nanmean(squeeze(nanmean(cleandat(:,:,:,2,2,6),3))-squeeze(nanmean(cleandat(:,:,:,1,2,6),3)));
fc = fc(idx);






























