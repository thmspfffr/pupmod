%% PLOT NUMBER OF ALTERED CORRELATIONS, INCLUDING STATISTICS
% --------------------------
% This script obtains the empircal number of altered correlations, by
% calling pupmod_compute_altered_correlations.m and obtains a corrected 
% p-values from a permutation distribution (computed in
% pupmod_src_powcorr_permtest.m). The actual p-values are obtained calling
% the function pupmod_all_powcorr_getstatistics.m)
% --------------------------
% CONTENTS
% --------------
% (1) PLOT: No stats
% (2) PLOT: P-Values (corrected) 
% (3) PLOT: Altered correlations
% (4) Circular graphs
% (5) Plot altered correlations (per voxel)
% --------------

clear

% -------------
% version of cleaned data: 
% v1: AAL, v2: 400 vertices (cortex)
% -------------
v = 1;
% -------------

% load  data
cd ~/pupmod/matlab/
cleandat = pupmod_loadpowcorr(v,1);

load redblue.mat
para = [];
para.n = 90;
para.transfer = 'to_bcn';

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

for ifoi = 1 : 13
  [h,~,~,s] = ttest(cleandat(:,:,:,1,2,ifoi),cleandat(:,:,:,1,1,ifoi),'dim',3);

  hh = tp_match_aal(para,h);
  s.tstat = tp_match_aal(para,s.tstat);
  
  hh = hh(include_bcn,include_bcn);
  s.tstat = s.tstat(include_bcn,include_bcn);
%   hh(hh==0)=-1;

  hh1 = triu((hh>0)&(s.tstat>0),1);
  hh2 = -tril((hh>0)&(s.tstat<0),1);

  h = hh1+hh2;
  figure; set(gcf,'color','w');
  imagesc(h,[-1 1]); axis off; colormap(redblue); axis square
  print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_aal_taskvsrest_f%d.eps',ifoi))

end