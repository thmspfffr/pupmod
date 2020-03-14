
clear 

SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

para          = [];
para.nfreq    = 1:17; % freqs 1 - 13
para.alpha    = 0.05; % alpha for adjacency
para.nperm    = 1000;
para.nsubs    = 25;
para.type     = 'global';
para.allperms = para.nperm/para.nsubs;
para.cleaned  = 0;
para.correction_method = 'single_threshold';
% ---------
% Obtain stats for atomoxetine condition
para.cond     = 'atx';

v = 1;
fc = pupmod_loadpowcorr(v,SUBJLIST,1); 
para1 = [];
para1.nfreq = 1 : 17;
para1.alpha = 0.05;
emp1 = pupmod_compute_altered_correlations(fc,para1);
emp1.reg = 0.05;
fc1 = squeeze(nanmean(nanmean(fc,1),2));

para.ver      = v;
para.emp      = emp1;
[outp_atx1]    = pupmod_src_powcorr_getstatistics(para);

v = 2;
fc = pupmod_loadpowcorr(v,SUBJLIST,1);
para1 = [];
para1.nfreq = 1 : 17;
para1.alpha = 0.05;
emp2 = pupmod_compute_altered_correlations(fc,para);
emp2.reg = 0.15;
fc2 = squeeze(nanmean(nanmean(fc,1),2));
para.ver      = v;
para.emp      = emp2;
[outp_atx2]    = pupmod_src_powcorr_getstatistics(para);

% v = 3;
% fc = pupmod_loadpowcorr(v,SUBJLIST,1); 
% 
% para1 = [];
% para1.nfreq = 1 : 17;
% para1.alpha = 0.05;
% emp3 = pupmod_compute_altered_correlations(fc,para);
% emp3.reg = 0.3;
% fc3 = squeeze(nanmean(nanmean(fc,1),2));
% para.ver      = v;
% para.emp      = emp3;
% [outp_atx3]    = pupmod_src_powcorr_getstatistics(para);

v = 4;
fc = pupmod_loadpowcorr(v,SUBJLIST,1); 

para1 = [];
para1.nfreq = 1 : 17;
para1.alpha = 0.05;
emp4 = pupmod_compute_altered_correlations(fc,para);
emp4.reg = 1;
fc4 = squeeze(nanmean(nanmean(fc,1),2));
para.ver      = v;
para.emp      = emp4;
[outp_atx4]    = pupmod_src_powcorr_getstatistics(para);

%%

foi_range = 2.^(2:.25:6);
figure_w;

subplot(4,2,1);
hold on
plot(emp1.n_n_atx(:,1),'b:')
plot(emp1.n_n_atx(:,2),'b-')
plot(emp1.n_p_atx(:,1),'r:')
plot(emp1.n_p_atx(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
ylabel('FAC');
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

title(sprintf('alpha0 = %.3f',emp1.reg))

subplot(4,2,2);
hold on
plot(emp1.n_n_dpz(:,1),'b:')
plot(emp1.n_n_dpz(:,2),'b-')
plot(emp1.n_p_dpz(:,1),'r:')
plot(emp1.n_p_dpz(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

subplot(4,2,3);
hold on
plot(emp2.n_n_atx(:,1),'b:')
plot(emp2.n_n_atx(:,2),'b-')
plot(emp2.n_p_atx(:,1),'r:')
plot(emp2.n_p_atx(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
ylabel('FAC');
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)
title(sprintf('alpha0 = %.3f',emp2.reg))

subplot(4,2,4);
hold on
plot(emp2.n_n_dpz(:,1),'b:')
plot(emp2.n_n_dpz(:,2),'b-')
plot(emp2.n_p_dpz(:,1),'r:')
plot(emp2.n_p_dpz(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

subplot(4,2,5);
hold on
plot(emp3.n_n_atx(:,1),'b:')
plot(emp3.n_n_atx(:,2),'b-')
plot(emp3.n_p_atx(:,1),'r:')
plot(emp3.n_p_atx(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
ylabel('FAC')
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp3.reg))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

subplot(4,2,6);
hold on
plot(emp3.n_n_dpz(:,1),'b:')
plot(emp3.n_n_dpz(:,2),'b-')
plot(emp3.n_p_dpz(:,1),'r:')
plot(emp3.n_p_dpz(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

subplot(4,2,7);
hold on
plot(emp4.n_n_atx(:,1),'b:')
plot(emp4.n_n_atx(:,2),'b-')
plot(emp4.n_p_atx(:,1),'r:')
plot(emp4.n_p_atx(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
ylabel('FAC'); xlabel('Frequency [Hz]'); 
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp4.reg))
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

subplot(4,2,8);
hold on
plot(emp4.n_n_dpz(:,1),'b:')
plot(emp4.n_n_dpz(:,2),'b-')
plot(emp4.n_p_dpz(:,1),'r:')
plot(emp4.n_p_dpz(:,2),'r-')

axis([0 18 0 0.3]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
xlabel('Frequency [Hz]'); 
set(gca,'ytick',0:0.1:0.3,'yticklabel',0:0.1:0.3)

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_compareregparameter_FAC.eps'))


[h11,p11]=ttest(squeeze(fc1(:,1,1,:)),squeeze(fc1(:,2,1,:)),'dim',1);
[h12,p12]=ttest(squeeze(fc1(:,1,2,:)),squeeze(fc1(:,2,2,:)),'dim',1);
[h21,p21]=ttest(squeeze(fc2(:,1,1,:)),squeeze(fc2(:,2,1,:)),'dim',1);
[h22,p22]=ttest(squeeze(fc2(:,1,2,:)),squeeze(fc2(:,2,2,:)),'dim',1);
[h31,p31]=ttest(squeeze(fc3(:,1,1,:)),squeeze(fc3(:,2,1,:)),'dim',1);
[h32,p32]=ttest(squeeze(fc3(:,1,2,:)),squeeze(fc3(:,2,2,:)),'dim',1);
[h41,p41]=ttest(squeeze(fc4(:,1,1,:)),squeeze(fc4(:,2,1,:)),'dim',1);
[h42,p42]=ttest(squeeze(fc4(:,1,2,:)),squeeze(fc4(:,2,2,:)),'dim',1);

figure_w;
subplot(4,2,1);
hold on
plot(squeeze(nanmean(fc1(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc1(:,2,1,:))),'k-')
plot(squeeze(nanmean(fc1(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc1(:,2,2,:))),'r-')

plot(find(p11<0.05),0.07*ones(sum(p11<0.05),1),'b*')
plot(find(p12<0.05),0.01*ones(sum(p12<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
ylabel('Mean corr.');
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp1.reg))

subplot(4,2,3);
hold on
plot(squeeze(nanmean(fc2(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc2(:,2,1,:))),'k-')
plot(squeeze(nanmean(fc2(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc2(:,2,2,:))),'r-')

plot(find(p21<0.05),0.07*ones(sum(p21<0.05),1),'b*')
plot(find(p22<0.05),0.01*ones(sum(p22<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
ylabel('Mean corr.');
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp2.reg))

subplot(4,2,5);
hold on
plot(squeeze(nanmean(fc3(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc3(:,2,1,:))),'k-')
plot(squeeze(nanmean(fc3(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc3(:,2,2,:))),'r-')

plot(find(p31<0.05),0.07*ones(sum(p31<0.05),1),'b*')
plot(find(p32<0.05),0.01*ones(sum(p32<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
ylabel('Mean corr.');
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp3.reg))

subplot(4,2,7);
hold on
plot(squeeze(nanmean(fc4(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc4(:,2,1,:))),'k-')
plot(squeeze(nanmean(fc4(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc4(:,2,2,:))),'r-')

plot(find(p41<0.05),0.07*ones(sum(p41<0.05),1),'b*')
plot(find(p42<0.05),0.01*ones(sum(p42<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
ylabel('Mean corr.'); xlabel('Frequency [Hz]'); 
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
title(sprintf('alpha0 = %.3f',emp4.reg))

[h11,p11]=ttest(squeeze(fc1(:,1,1,:)),squeeze(fc1(:,3,1,:)),'dim',1);
[h12,p12]=ttest(squeeze(fc1(:,1,2,:)),squeeze(fc1(:,3,2,:)),'dim',1);
[h21,p21]=ttest(squeeze(fc2(:,1,1,:)),squeeze(fc2(:,3,1,:)),'dim',1);
[h22,p22]=ttest(squeeze(fc2(:,1,2,:)),squeeze(fc2(:,3,2,:)),'dim',1);
[h31,p31]=ttest(squeeze(fc3(:,1,1,:)),squeeze(fc3(:,3,1,:)),'dim',1);
[h32,p32]=ttest(squeeze(fc3(:,1,2,:)),squeeze(fc3(:,3,2,:)),'dim',1);
[h41,p41]=ttest(squeeze(fc4(:,1,1,:)),squeeze(fc4(:,3,1,:)),'dim',1);
[h42,p42]=ttest(squeeze(fc4(:,1,2,:)),squeeze(fc4(:,3,2,:)),'dim',1);

subplot(4,2,2);

hold on
plot(squeeze(nanmean(fc1(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc1(:,3,1,:))),'k-')
plot(squeeze(nanmean(fc1(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc1(:,3,2,:))),'r-')

plot(find(p11<0.05),0.07*ones(sum(p11<0.05),1),'k*')
plot(find(p12<0.05),0.01*ones(sum(p12<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))

subplot(4,2,4);
hold on
plot(squeeze(nanmean(fc2(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc2(:,3,1,:))),'k-')
plot(squeeze(nanmean(fc2(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc2(:,3,2,:))),'r-')

plot(find(p21<0.05),0.07*ones(sum(p21<0.05),1),'k*')
plot(find(p22<0.05),0.01*ones(sum(p22<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))

subplot(4,2,6);
hold on
plot(squeeze(nanmean(fc3(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc3(:,3,1,:))),'k-')
plot(squeeze(nanmean(fc3(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc3(:,3,2,:))),'r-')

plot(find(p31<0.05),0.07*ones(sum(p31<0.05),1),'k*')
plot(find(p32<0.05),0.01*ones(sum(p32<0.05),1),'r*')


axis([0 18 0 0.075]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))

subplot(4,2,8);
hold on
plot(squeeze(nanmean(fc4(:,1,1,:))),'k:')
plot(squeeze(nanmean(fc4(:,3,1,:))),'k-')
plot(squeeze(nanmean(fc4(:,1,2,:))),'r:')
plot(squeeze(nanmean(fc4(:,3,2,:))),'r-')

plot(find(p41<0.05),0.07*ones(sum(p41<0.05),1),'k*')
plot(find(p42<0.05),0.01*ones(sum(p42<0.05),1),'r*')

axis([0 18 0 0.075]); tp_editplots
set(gca,'xtick',1:4:17,'xticklabel',foi_range(1:4:17))
xlabel('Frequency [Hz]'); 

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_powcorr_compareregparameter.eps'))

%% SOURCE LEVEL PLOTS

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

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

close all

for v = [1 2 4]
  
  ifoi = [6 7 8 9];
  icond = 2;

  cmap  = autumn;
  eval(sprintf('par   = nanmean(emp%d.n_p_atx_pervoxel(:,ifoi,icond),2);',v));
  prctile(par,95)
  par(par<prctile(par,95))=0;
  
%   par(outp_atx2.pval_vox_p_atx(:,icond,ifoi)>0.05)=0;
  
  cmap      = [cmap; 0.98*ones(1,3); cmap];
  para      = [];
  para.clim = [-max(par) max(par)]
  para.cmap = cmap;
  para.grid = grid;
  para.dd   = 0.75;
%   para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_atx_compareregparameter_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
%   tp_plot_surface(par,para)
  
end
%%
% for v = 1:4
%   
%   ifoi = [12:14];
%   icond = 2;
% 
%   cmap  = autumn;
%   eval(sprintf('par   = nanmean(emp%d.n_p_atx_pervoxel(:,ifoi,icond),2);',v));
%   prctile(par,95)
%   par(par<prctile(par,95))=0;
%   cmap      = [cmap; 0.98*ones(1,3); cmap];
%   para      = [];
%   para.clim = [-max(par) max(par)]
%   para.cmap = cmap;
%   para.grid = grid;
%   para.dd   = 0.75;
% %   para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_atx_compareregparameter_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
% %   tp_plot_surface(par,para)
% %   
% end


for v = [1 2 4]
  
  ifoi = [6 7 8 9 10];
  icond = 1;

  cmap      = autumn;
  cmap(:,1) = 0; 
  cmap(:,3) = 1;
  eval(sprintf('par   = nanmean(emp%d.n_n_dpz_pervoxel(:,ifoi,icond),2);',v));
  prctile(par,95)
%   par(outp_atx2.pval_vox_n_dpz(:,icond,ifoi)>0.025)=0;
    par(par<prctile(par,95))=0;

  cmap      = [cmap; 0.98*ones(1,3); cmap];
  para      = [];
  para.clim = [-max(par) max(par)]
  para.cmap = cmap;
  para.grid = grid;
  para.dd   = 0.75;
%   para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_dpz_compareregparameter_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
%   tp_plot_surface(par,para)
  
end

% 
% for v = 1:4
%   
%   ifoi = [1:2];
%   icond = 1;
% 
%   cmap      = autumn;
%   cmap(:,1) = 0; 
%   cmap(:,3) = 1;
%   eval(sprintf('par   = nanmean(emp%d.n_n_dpz_pervoxel(:,ifoi,icond),2);',v));
%     prctile(par,95)
% 
%   par(par<prctile(par,95))=0;
%   cmap      = [cmap; 0.98*ones(1,3); cmap];
%   para      = [];
%   para.clim = [-max(par) max(par)]
%   para.cmap = cmap;
%   para.grid = grid;
%   para.dd   = 0.75;
% %   para.fn   = sprintf('~/pupmod/plots/pupmod_plot_alteredcorr_dpz_compareregparameter_f%s_c%d_v%d.png',regexprep(num2str(ifoi),' ',''),icond,v);
% %   tp_plot_surface(par,para)
%   
% end

%% Correlate effects

for isubj = 1 : 28
  for ifoi = 1 : 17
  tmp1 = squeeze(fc(:,:,isubj,3,1,ifoi)) - squeeze(fc(:,:,isubj,1,1,ifoi));
  tmp2 = squeeze(fc(:,:,isubj,2,2,ifoi)) - squeeze(fc(:,:,isubj,1,2,ifoi));
  r(isubj,ifoi)=corr(tmp1(mask),tmp2(mask))
  end
end



%%
fc = pupmod_loadpowcorr(1,SUBJLIST,0);

high_in_block1 = fc(:,:,:,1,2,:,1)>nanmean(nanmean(fc(:,:,:,1,2,:,1),1),2);
fc2 = fc(:,:,:,2,2,:,2)-fc(:,:,:,1,2,:,2);
mean(fc2(high_in_block1))
nanmean(fc2(~high_in_block1))