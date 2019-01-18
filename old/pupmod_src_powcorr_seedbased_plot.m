v = 12;

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v))

restoredefaultpath

addpath ~/pconn/matlab
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
% addpath /home/gnolte/neuconn/OLD/matlab/rest/
addpath(genpath('/home/gnolte/meth'));

%%

[h,p,~,s] = ttest(cleandat(:,:,:,2,1,7),cleandat(:,:,:,1,1,7),'dim',3);

for i = 1 : 90 
  for j = 1 : 90
    
    p(i,j) = 1-sum(abs(s.tstat(i,j))>max_cnt1_atx)/nperm;
    
  end
end
    

%% SEED BASED CONNECTIVITY
figure;
load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',4,1,1,9))
grid = sa.grid_cortex_lowres;

% load('~/pupmod/proc/grid_cortexlow.mat')

% VISUAL
%   vis = [-47 -69 -3; 47 -69 3]./10;
% plt = 'vis';
% 
% SOMATOSENSORY
% vis = [-42 -26 54; 38 -32 48]./10;
% plt = 'som';

% AUDITORY
vis = [-54 -22 10; 52 -24 12]./10;
plt = 'aud';

% vis = [-39 -54 32; 52 -24 12]./10;
% plt = 'lpc';

% vis = [-3 39 -2; 3 39 -2]./10;
% plt = 'mpfc';

plot3(grid(:,1),grid(:,2),grid(:,3),'y.'); hold on

for i = 1 : length(grid)
  d1(i) = abs(sqrt([grid(i,1)-vis(1,1)].^2+[grid(i,2)-vis(1,2)].^2+[grid(i,3)-vis(1,3)].^2));
  d2(i) = abs(sqrt([grid(i,1)-vis(2,1)].^2+[grid(i,2)-vis(2,2)].^2+[grid(i,3)-vis(2,3)].^2));
end
[~,i1]=min(d1)

plot3(grid(i1,1),grid(i1,2),grid(i1,3),'k.','markersize',50)
plot3(vis(1,1),vis(1,2),vis(1,3),'r.','markersize',50)
plot3(vis(2,1),vis(2,2),vis(2,3),'r.','markersize',50)
[~,i2]=min(d2)

plot3(grid(i2,1),grid(i2,2),grid(i2,3),'b.','markersize',50)
load sa_meg_template;
%%



para.plot = 'i1';
ifoi = 8;


if strcmp(para.plot,'i2')
  ii = i2;
else
  ii = i1;
end

pp = para.plot;

clear h p

addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/

avg_corr = squeeze(nanmean(nanmean(cleandat(ii,:,:,1,1,ifoi),4),2));

for i = 1 : 400
  [h(i),p(i)] = ttest(squeeze(mean(cleandat(ii,i,:,1,1,ifoi),4)),avg_corr);
end

sign = p<.05;
para = [];


m = squeeze(mean(mean(cleandat(ii,:,:,1,1,ifoi),4),3));
sign(isnan(m))=[];
m(isnan(m)) = [];

grid = sa.grid_cortex_lowres;

% roi = [grid(i1,:);grid(i2,:)] ;
roi = [grid(ii,:); [-grid(ii,1) grid(ii,2) grid(ii,3)]] ;

grid(i1,:) = [];

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = .5;
% m2 = spatfiltergauss(m,g1,dd,g2);

z2 = spatfiltergauss(m(:),g1,dd,g2).*(1/0.577);
% para.sign = spatfiltergauss(sign(:),g1,dd,g2);

para = [] ;
para.colorlimits = [0.080 0.10];

viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
cmap = inferno;
cmap(1:100,:)=0.98*ones(100,3);
para.mymarkersize = 40;
% para.fsv = 10
tp_showsource(z2,cmap,sa_meg_template,para);

print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_src_powcorr_seeds_cortex_%s_%s_f%d_v%d.jpg',plt,pp,ifoi,v))

% SHOW MARKER -------------------------------------------------------------
grid = sa.grid_cortex_lowres;
roi = [grid(ii,:); [-grid(ii,1) grid(ii,2) grid(ii,3)]] ;
m = zeros(400,1);

m(i1,:)  = 1;
m(i2,:)  = 1;

dd          = 0.1;
z2          = spatfiltergauss(m(:),grid,dd,g2);

para = [] ;
para.colorlimits = [0.99 1];

viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
cmap = inferno;
% cmap(1:100,:)=0.98*ones(100,3);
tp_showsource(z2,cmap,sa_meg_template,para);

print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_src_powcorr_seeds_cortex_marker_%s_%s_f%d_v%d.jpg',plt,pp,ifoi,v))
