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
for icond = 1 : 2
  for ifoi = 1:13
%     
    if ~exist(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d_processing.txt'],icond,ifoi,vv))
      system(['touch ' outdir sprintf('pupmod_all_src_degree_c%d_f%d_v%d_processing.txt',icond,ifoi,vv)]);
    else
      continue
    end
%     
    j = 1:size(cleandat,1);

    for il = 1 : size(cleandat,1)

      jj = j(j~=il);

      for jl = 1 : size(cleandat,1)

        jjj = j(j~=jl);

        fc_tmp = cleandat(:,:,:,1:3,icond,ifoi);

        fprintf('Computing location %d %d ...\n',il,jl);

        x = squeeze(abs(fc_tmp(il,jl,:,:)));
        
        % Reference connections:
        % il -> all others
        % jl -> all others
        % Reference is currently determined from placebo condition       
        x_ref1 = abs(squeeze(nanmean(fc_tmp(il,jj,:,1),4)));
        x_ref2 = abs(squeeze(nanmean(fc_tmp(jl,jjj,:,1),4)));

        % Compute mean and std over reference connections
        m1 = mean(x_ref1); s1 = std(x_ref1);
        m2 = mean(x_ref2); s2 = std(x_ref2);

        % Compare connection to reference, obtain z-score
        z11 = (abs(x(:,1))' - m1)./s1;
        z12 = (abs(x(:,1))' - m2)./s2;

        z21 = (abs(x(:,2))' - m1)./s1;
        z22 = (abs(x(:,2))' - m2)./s2;

        z31 = (abs(x(:,3))' - m1)./s1;
        z32 = (abs(x(:,3))' - m2)./s2;
        
        % Test against zero
        % Left tailed: count only those correlations that
        % are significantly LARGER than zero.
        
        z = [z11; z12; z21; z22; z31; z32];
        [~,p] = ttest(zeros(size(z)),z,'tail','left','dim',2);

        th11(jl,:) = p(1) < 0.01/2;
        th12(jl,:) = p(2) < 0.01/2;
        th21(jl,:) = p(3) < 0.01/2;
        th22(jl,:) = p(4) < 0.01/2;
        th31(jl,:) = p(5) < 0.01/2;
        th32(jl,:) = p(6) < 0.01/2;

      end

      % any connection significant?
      th1(il,jj,:) = (th11(jj,:) + th12(jj,:)) > 0; clear th11 th12
      th2(il,jj,:) = (th21(jj,:) + th22(jj,:)) > 0; clear th21 th22
      th3(il,jj,:) = (th31(jj,:) + th32(jj,:)) > 0; clear th31 th32

    end

    k(:,1) = sum(th1)/size(cleandat,1);
    k(:,2) = sum(th2)/size(cleandat,1);
    k(:,3) = sum(th3)/size(cleandat,1);

    save(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,vv),'k');
  
    clear k th1 th2 th3 th11 th12 th21 th22 th31 th32
    
  end
end

error('!')

% Add more metrics

%% PLOT RESULTS OF DEGREE COMPUTATION
foi_range = unique(round(2.^[1:.5:7]))

figure; hold on

for icond = 1 : 2

v = 121

for ifoi = 1 : 13
  
  load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));
  
  kk(:,:,ifoi) = k; clear k
  
end
  


subplot(2,1,icond); hold on

plot(log10(foi_range),squeeze(mean(kk(:,1,:),1)).*100,'color',[.7 .7 .7],'linewidth',2)
plot(log10(foi_range),squeeze(mean(kk(:,2,:),1)).*100,'color',[1 0.2 0.2],'linewidth',2)
plot(log10(foi_range),squeeze(mean(kk(:,3,:),1)).*100,'color',[0.2 0.2 1],'linewidth',2)

xlabel('Frequency [Hz]'); ylabel('Degree [%]')

set(gca,'tickdir','out','xtick',log10(foi_range([1 3 5 7 9 11 13])),'xticklabel',[2 4 8 16 32 64 128])
% 
end
  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_degree_spectrum.pdf'))

  figure;
  
for icond = 1 : 2

  v = 121

  for ifoi = 1 : 13

    load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));

    kk(:,:,ifoi) = k; clear k

  end
  
  subplot(2,1,icond)
  
  deg = mean(squeeze(mean(kk(:,:,:),1)),2);
  bar(deg)
  
  axis([0 4 0 0.15]); ylabel('Mean degree [%]')
  axis square; %tp_editplots
  
end

  print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_mean_degree.pdf'))

  
%% PLOT ON SURFACE
v = 121
% ifoi = 7;
% icond = 2;
  for icond = 1

for ifoi = 5
load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));

var2plot = k(:,1);

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
  
  
  
  

  
  

