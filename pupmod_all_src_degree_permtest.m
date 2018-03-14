%% pupmod_all_src_degree
% plot degree of cleaned signal
% obtain cleanined signal from pupmod_all_src_peripheral*** *check)
% Goal: replicate hipp et al., nn

clear

v = 12;

outdir = '~/pupmod/proc/conn/';
  
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));


%%

s_res = squeeze(cleandat(:,:,:,:,1,:));
s_cnt = squeeze(cleandat(:,:,:,:,2,:));

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')

nperm = 10000; alp = 0.05;

par.subs = 50;
par.allperms = nperm/par.subs;

if ~exist(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v),'all_idx1');
else
  load(sprintf('~/pupmod/proc/pupmod_src_degree_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v));
end

dat_cnt1 = squeeze(cleandat(:,:,:,[1 2],2,:)); 
dat_res1 = squeeze(cleandat(:,:,:,[1 2],1,:)); 
dat_cnt2 = squeeze(cleandat(:,:,:,[1 3],2,:));
dat_res2 = squeeze(cleandat(:,:,:,[1 3],1,:));

for iperm = 1 : par.allperms
  
  if ~exist(sprintf(['~/pupmod/proc/' 'pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' '~/pupmod/proc/' sprintf('pupmod_src_degree_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    

for icond = 1 : 2
  for ifoi = 1:13
    
    if ~exist(sprintf([outdir 'pupmod_all_src_degree_permtest_c%d_f%d_v%d_processing.txt'],icond,ifoi,v))
      system(['touch ' outdir sprintf('pupmod_all_src_permtest_degree_c%d_f%d_v%d_processing.txt',icond,ifoi,v)]);
    else
      continue
    end
    
    j = 1:size(cleandat,1);

    for il = 1 : size(cleandat,1)

      jj = j(j~=il);

      for jl = 1 : size(cleandat,1)

        jjj = j(j~=jl);

        fc_tmp = cleandat(:,:,:,1:3,icond,ifoi);

        fprintf('Computing location %d %d ...\n',il,jl);

        x = squeeze(fc_tmp(il,jl,:,:));

        x_ref1 = squeeze(nanmean(fc_tmp(il,jj,:,:),4));
        x_ref2 = squeeze(nanmean(fc_tmp(jjj,jl,:,:),4));

        m1 = mean(x_ref1); s1 = std(x_ref1);
        m2 = mean(x_ref2); s2 = std(x_ref2);

        z11 = (x(:,1)' - m1)./s1;
        z12 = (x(:,1)' - m2)./s2;

        z21 = (x(:,2)' - m1)./s1;
        z22 = (x(:,2)' - m2)./s2;

        z31 = (x(:,3)' - m1)./s1;
        z32 = (x(:,3)' - m2)./s2;

        [~,p11] = ttest(z11);
        [~,p12] = ttest(z12);

        [~,p21] = ttest(z21);
        [~,p22] = ttest(z22);

        [~,p31] = ttest(z31);
        [~,p32] = ttest(z32);

        th11(jl,:) = p11 < 0.05/2;
        th12(jl,:) = p12 < 0.05/2;
        th21(jl,:) = p21 < 0.05/2;
        th22(jl,:) = p22 < 0.05/2;
        th31(jl,:) = p31 < 0.05/2;
        th32(jl,:) = p32 < 0.05/2;

      end

      th11 = th11(jj,:);
      th12 = th12(jj,:);
      th21 = th21(jj,:);
      th22 = th22(jj,:);
      th31 = th31(jj,:);
      th32 = th32(jj,:);

      % any connection significant?
      th1(il,jj,:) = (th11 + th12) > 0; clear th11 th12
      th2(il,jj,:) = (th21 + th22) > 0; clear th21 th22
      th3(il,jj,:) = (th31 + th32) > 0; clear th31 th32

    end

    k(:,1) = sum(th1)/size(cleandat,1);
    k(:,2) = sum(th2)/size(cleandat,1);
    k(:,3) = sum(th3)/size(cleandat,1);

    save(sprintf([outdir 'pupmod_all_src_degree_permtest_c%d_f%d_v%d.mat'],icond,ifoi,v),'k');
  
  end
end

error('!')

% Add more metrics

%% PLOT RESULTS OF DEGREE COMPUTATION
foi_range = unique(round(2.^[1:.5:7]))

icond = 1;

for ifoi = 1 : 13
  
  load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));
  
  kk(:,:,ifoi) = k; clear k
  
end
  


figure; hold on
% subplot(1,2,1); 

plot(log10(foi_range),squeeze(mean(kk(:,1,:),1)).*100,'color',[.7 .7 .7],'linewidth',2)
plot(log10(foi_range),squeeze(mean(kk(:,2,:),1)).*100,'color',[1 0.5 0.2],'linewidth',2)
plot(log10(foi_range),squeeze(mean(kk(:,3,:),1)).*100,'color',[0.2 0.5 1],'linewidth',2)

xlabel('Frequency [Hz]'); ylabel('Degree [%]')

set(gca,'tickdir','out','xtick',log10(foi_range([1 3 5 7 9 11 13])),'xticklabel',[2 4 8 16 32 64 128])
% 
% box off; axis square; axis([1 14 0 35])
% 
% subplot(1,2,2); hold on
% 
% plot(k(1,:,2)'.*100,'color',[.7 .7 .7],'linewidth',2)
% plot(k(2,:,2)'.*100,'color',[1 0.5 0.2],'linewidth',2)
% plot(k(3,:,2)'.*100,'color',[0.2 0.5 1],'linewidth',2)
% 
% xlabel('Frequency [Hz]'); ylabel('Degree [%]')
% 
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% 
% box off; axis square; axis([1 14 0 35])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_degree.pdf'))

%% PLOT ON SURFACE

ifoi = 6;
icond = 2;
load(sprintf([outdir 'pupmod_all_src_degree_c%d_f%d_v%d.mat'],icond,ifoi,v));

var2plot = k(:,3)-k(:,1);

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

para = [] ;
para.colorlimits = [min(var2plot) max(var2plot)];

viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
cmap = jet;
% cmap(1:100,:)=0.98*ones(100,3);

tp_showsource(z2,cmap,sa_meg_template,para);



% 
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
  
  
  
  

  
  

