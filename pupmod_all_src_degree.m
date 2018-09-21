%% pupmod_all_src_degree
% plot degree of cleaned signal
% obtain cleanined signal from pupmod_all_src_peripheral*** *check)
% Goal: replicate hipp et al., nn

clear

vv = 121;
v = 12;

outdir = '~/pupmod/proc/conn/';
  
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));



%%

fcsize = size(cleandat,1);
para = [];
para.alpha = 0.01;
para.nfreq = 13;


deg_atx = tp_degree(cleandat(:,:,:,[1 2],[1:2],:),para);
deg_dpz = tp_degree(cleandat(:,:,:,[1 3],[1:2],:),para);

deg_atx = squeeze(nansum(reshape(deg_atx,[fcsize^2 13 2]))/fcsize^2);
deg_dpz = squeeze(nansum(reshape(deg_dpz,[fcsize^2 13 2]))/fcsize^2);

deg_atx_task = tp_degree(cleandat(:,:,:,[1 2],[2],:),para);
deg_atx_task = squeeze(nansum(reshape(deg_atx_task,[fcsize^2 13 2]))/fcsize^2);

deg_dpz_task = tp_degree(cleandat(:,:,:,[1 3],[2],:),para);
deg_dpz_task = squeeze(nansum(reshape(deg_dpz_task,[fcsize^2 13 2]))/fcsize^2);


%% CONCATENATE PERMUTATONS AND COMPUTE STATS
v = 12;
nperm = 50000; 
par.subs = 100;
par.allperms = nperm/par.subs;

for iperm = 1 : par.allperms
  
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  
  perm_k_atx(:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_atx;
  perm_k_dpz(:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_dpz;
%   perm_k_atx_pervoxel(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_atx_pervoxel;
%   perm_k_dpz_pervoxel(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_dpz_pervoxel;
%   perm_k_tvr(:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_tvr;
%   perm_k_tvr_pervoxel(:,:,:,(iperm-1)*par.subs+1:(iperm)*par.subs) = outp.perm_k_tvr_pervoxel;
  
end

%% P-VALUES

for ifoi = 1 : 13
  p_atx_rest(ifoi) = 1-sum(deg_atx(ifoi,2) > squeeze(perm_k_atx(ifoi,2,:)))/50000;
  p_dpz_task(ifoi) = 1-sum(deg_dpz(ifoi,1) > squeeze(perm_k_dpz(ifoi,1,:)))/50000;
end

figure; hold on
plot(p_atx_rest)
plot(p_atx_task)
plot(p_dpz_rest)
plot(p_dpz_task)


%% PLOT ON SURFACE
v = 121
% ifoi = 7;
% icond = 2;
  for icond = 1

for ifoi = 5
load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));

var2plot = k(:,1);

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/sa_template.mat;
  sa_meg_template = sa_template;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
end

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = .5;

z2 = spatfiltergauss(var2plot(:),g1,dd,g2);
% 
% para = [] ;
para.colorlimits = [min(var2plot) max(var2plot)];
% % para.colorlimits = [-max(abs(var2plot)) max(abs(var2plot))];
% 
% viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
% % 
igr = 80;
cmap = inferno;
c = [linspace(0.98,cmap(igr,1),igr)' linspace(0.98,cmap(igr,2),igr)' linspace(0.98,cmap(igr,3),igr)'];
cmap(1:igr,:)=c;
para.filename = sprintf('~/pupmod/plots/pupmod_src_degree_f%d_c%d_v%d.png',ifoi,icond,v)
% 
tp_showsource(z2,cmap,sa_meg_template,para);

  end
  end

  %% PLOT CONTRAST
 
  igr = 50;
  
% cmap1 = cbrewer('seq','YlOrRd',126); cmap1 = cmap1(1:100,:); cmap1 = cmap1(end:-1:1,:);
% ccmap1 = [linspace(0.98,cmap1(1,1),igr)' linspace(0.98,cmap1(1,2),igr)' linspace(0.98,cmap1(1,3),igr)'];
% ccmap1 = ccmap1(1:end-1,:); ccmap1 = ccmap1(end:-1:1,:);
% cmap1 = [ccmap1(end:-1:1,:); cmap1];
% 
% cmap2 = cbrewer('seq','YlGnBu',126); cmap2 = cmap2(1:100,:); cmap2 = cmap2(end:-1:1,:);
% ccmap2 = [linspace(0.98,cmap2(1,1),igr)' linspace(0.98,cmap2(1,2),igr)' linspace(0.98,cmap2(1,3),igr)'];
% ccmap2 = ccmap2(1:end-1,:); ccmap2 = ccmap2(end:-1:1,:);
% cmap2 = [ccmap2(end:-1:1,:); cmap2];
% cmap = [cmap2(end:-1:1,:); ones(20,3); cmap1]

clear cmap1 cmap2 cmap autumn
% 
cmap1 = hot(150);
cmap2 = cmap1(:,[3 2 1]);
cmap = [cmap2(end-20:-1:30,:); .98*ones(30,3);  cmap1(30:end-20,:)];

% igr = 80;
% cmap = inferno;
% c = [linspace(0.98,cmap(igr,1),igr)' linspace(0.98,cmap(igr,2),igr)' linspace(0.98,cmap(igr,3),igr)'];
% cmap(1:igr,:)=c;
% 
% cmap = [cmap(end:-1:1,[3 2 1]); cmap]

%%

v = 121

foi = [1:13];
cond = [1:2];
icontr = 2;

for ifoi = foi
  for icond = cond

load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));

if icontr == 1
  var2plot = k(:,2) -k(:,1);
elseif icontr == 2
  var2plot = k(:,3) -k(:,1);
else
  var2plot = k(:,2) -k(:,3);
end

if ~exist('sa_meg_template','var')
  load sa_meg_template;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
end

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = .5;

z2 = spatfiltergauss(var2plot(:),g1,dd,g2);
% 
para = [] ;
% para.colorlimits = [min(var2plot) max(var2plot)];
% para.colorlimits = [-max(abs(var2plot)) max(abs(var2plot))];
para.colorlimits = [-0.45 0.45];

para.viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
% 
% igr = 80;
% cmap = cbrewer('div', 'RdBu', 256,'pchip');
% cmap = cmap(end:-1:1,:);

para.filename = sprintf('~/pupmod/plots/pupmod_src_degree_contr%d_f%d_c%d_v%d.png',icontr,ifoi,icond,v)

tp_showsource(z2,cmap,sa_meg_template,para);

  end
end

% for icond = 1 : 2
% for ifoi = 1:13
% j = 1:90;
% 
% 
% for il = 1 : size(s_fc,1)
%   
%   jj = j(j~=il);
%   
%   for jl = 1 : size(s_fc,1)
%     
%     jjj = j(j~=jl);
%     
%     fc_tmp = s_fc(:,:,:,1:3,icond,ifoi);
%     
%     fprintf('Computing location %d %d ...\n',il,jl);
%     
%     x = squeeze(fc_tmp(il,jl,:,:));
%     
%     x_ref1 = squeeze(nanmean(fc_tmp(il,jj,:,:),4));
%     x_ref2 = squeeze(nanmean(fc_tmp(jjj,jl,:,:),4));
%     
%     m1 = mean(x_ref1); s1 = std(x_ref1);
%     m2 = mean(x_ref2); s2 = std(x_ref2);
%     
%     z11 = (x(:,1)' - m1)./s1;
%     z12 = (x(:,1)' - m2)./s2;
%     
%     z21 = (x(:,2)' - m1)./s1;
%     z22 = (x(:,2)' - m2)./s2;
%     
%     z31 = (x(:,3)' - m1)./s1;
%     z32 = (x(:,3)' - m2)./s2;
%        
%     [~,p11] = ttest(z11);
%     [~,p12] = ttest(z12);
%     
%     [~,p21] = ttest(z21);
%     [~,p22] = ttest(z22);
%     
%     [~,p31] = ttest(z31);
%     [~,p32] = ttest(z32);
%     
%     th11(jl,:) = p11 < 0.01/2;
%     th12(jl,:) = p12 < 0.01/2;
%     th21(jl,:) = p21 < 0.01/2;
%     th22(jl,:) = p22 < 0.01/2;
%     th31(jl,:) = p31 < 0.01/2;
%     th32(jl,:) = p32 < 0.01/2;
%     
%   end
%   
%   th11 = th11(jj,:);
%   th12 = th12(jj,:);
%   th21 = th21(jj,:);
%   th22 = th22(jj,:);
%   th31 = th31(jj,:);
%   th32 = th32(jj,:);
%   
%   % any connection significant?
%   th1(il,jj,:) = (th11 + th12) > 0; clear th11 th12
%   th2(il,jj,:) = (th21 + th22) > 0; clear th21 th22
%   th3(il,jj,:) = (th31 + th32) > 0; clear th31 th32
%     
% end
% 
% k(1,ifoi,icond) = sum(sum(th1))/8010
% k(2,ifoi,icond) = sum(sum(th2))/8010
% k(3,ifoi,icond) = sum(sum(th3))/8010
% 
% end
% 
% end

% for ifoi = 1 : 13
%   tmp = squeeze(s_fc(:,:,:,1,:,ifoi));
%   
%   [h,~,~,s] = ttest(tmp(:,:,:,2),tmp(:,:,:,1),'dim',3);
%   
%   n_p(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)>0))/4005;
%   n_n(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)<0))/4005;
%   
% end



% for iperm = 1 : nperm
  
  
  
  

  
  

