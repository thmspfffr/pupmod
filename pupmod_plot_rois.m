%%  PLOT ROIs on cortical surface

% DIMENSION 1 (in grid): left (negative) to right (positive)
% DIMENSION 2 (in grid): front (positive) to back (negative)
% DIMENSION 3 (in grid): top (positive) to bottom (negative)

% fc = pupmod_loadpowcorr(1,SUBJLIST,1);
% fc_rest = squeeze(fc(:,:,:,1,1,:));

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
istr = 4;

if istr == 1
  % 1 -----------------------
  str = 'supramarginal';
  r_roi = [[49 -31 50]./10; [60 -50 34]./10; [60 -31 34]./10];
  l_roi = [[-46 -35 56]./10; [-60 -50 30]./10; [-60 -35 30]./10];
  gridd{1} = [grid; r_roi(1,:); l_roi(1,:)];
  gridd{2} = [grid; r_roi(2,:); l_roi(2,:)];
  gridd{3} = [grid; r_roi(3,:); l_roi(3,:)];
  gridd{4} = gridd{3};
elseif istr == 2
  % 2 -----------------------
  str = 'intraparietal_sulc';
  r_roi = [[33 -53 60]./10; [33 -70 47]./10; [50 -53 47]./10];
  l_roi = [[-31 -54 55]./10; [-31 -70 46]./10;[-50 -54 46]./10];
  gridd{1} = [grid; r_roi(1,:); l_roi(1,:)];
  gridd{2} = [grid; r_roi(2,:); l_roi(2,:)];
  gridd{3} = [grid; r_roi(3,:); l_roi(3,:)];
  gridd{4} = gridd{3};
elseif istr == 3
  % 3 -----------------------
  str = 'v5';
  r_roi = [[45 -68 40]./10; [45 -78 1]./10; [52 -68 1]./10];
  l_roi = [[-47 -72 40]./10; [-47 -78 3]./10; [-52 -72 3]./10];
  gridd{1} = [grid; r_roi(1,:); l_roi(1,:)];
  gridd{2} = [grid; r_roi(2,:); l_roi(2,:)];
  gridd{3} = [grid; r_roi(3,:); l_roi(3,:)];
  gridd{4} = gridd{3};
elseif istr == 4
  % 4 -----------------------
  str = 'lpc';
  r_roi = [[46 -45 52]./10; [46 -60 39]./10; [60 -45 39]./10];
  l_roi = [[-39 -54 52]./10; [-39 -74 32]./10; [-59 -54 32]./10];
  gridd{1} = [grid; r_roi(1,:); l_roi(1,:)];
  gridd{2} = [grid; r_roi(2,:); l_roi(2,:)];
  gridd{3} = [grid; r_roi(3,:); l_roi(3,:)];
  gridd{4} = gridd{3};
elseif istr == 5
elseif istr == 6
elseif istr == 7
elseif istr == 8
elseif istr == 9
elseif istr == 10
elseif istr == 11
elseif istr == 12
  
end

par = zeros(400,1);
par(401) = 1;
par(402) = 1;

para = [];
para.clim   = [0 1];
para.cmap   = jet;
para.grid 	= gridd;
para.dd     = 0.5;
para.r_roi  = r_roi;
para.l_roi  = l_roi;
para.fn = sprintf('~/pupmod/plots/rois_%s.png',str)

tp_plot_roi(par,sa_template,para)

%% PLOT CORRELATION FROM  SEED TO THE REST

% SOM = [-42 -26 54; 38 -32 48]./10
SOM = [-20 -86 18; 16 -80 26]./10
% SOM = [-54 -22 10; 52 -24 12]./10


for i= 1 :400
  
  d_l(i) = sqrt((grid(i,1)-SOM(1,1))^2 + (grid(i,2)-SOM(1,2))^2 + (grid(i,3)-SOM(1,3))^2);
  d_r(i) = sqrt((grid(i,1)-SOM(2,1))^2 + (grid(i,2)-SOM(2,2))^2 + (grid(i,3)-SOM(2,3))^2);
  
end
[~,idx_l]=min(d_l);
[~,idx_r]=min(d_r);

idx = 1 : 400;
for i = 1 : 400
  
  others = squeeze(mean(fc(idx_l,~ismember(idx, i),:,9),2));
  this = squeeze(fc(idx_l,i,:,9));
  h(i)=ttest(this,others,'tail','right');
end


ifoi = 9; icond = 1;

cmap      = plasma;
para      = [];
para.clim = [0.05 0.10];
par = squeeze(nanmean((fc_rest(idx_l,:,:,ifoi)),3));
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn   = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)

%% MEAN CORRELATION



ifoi = 9; icond = 1;

cmap      = plasma;
para      = [];
para.clim = [.35 .43];
par = squeeze(nanmean(nanmean(fc(:,:,:,ifoi),2),3));
para.cmap = cmap;
para.grid = grid;
para.dd   = 0.75;
para.fn   = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)
