%% pupmod_src_powcorr_plot

clear all

% --------------------------------------------------------
% VERSION 10
% --------------------------------------------------------
v               = 9;
v_postproc      = 6;
v_grid          = 4; % 4 = aal
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal';
foi_range       = unique(round(2.^[1:.5:7]));
para.smo        = foi_range./4;
para.segleng    = 1 ./ para.smo;
lpc             = 0;
timevariant     = 1;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/';
addpath /home/tpfeffer/pconn/matlab/


%%
clear s s1 s2 fc_mean
v = 10;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

addpath ~/pconn/matlab/
  
ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
%         for iep = 1 : size(powcorr,3)
%           for jep = 1 : size(powcorr,3)
%             
%             fcd(iep,jep) = corr(nonzeros(triu(powcorr(:,:,iep),1)),nonzeros(triu(powcorr(:,:,jep),1)));
%           
%           end
%         end
%         
%         [h(1,:,m,ifoi,isubj,iblock),n]=hist(nonzeros(triu(fcd,1)),-1:0.01:1);

  s(:,:,isubj,m,ifoi,iblock) =  powcorr;
              
      end
    end
  end
end

save('~/pupmod/proc/pupmod_src_fcd.mat','h')

s = squeeze(nanmean(s(:,:,SUBJLIST,:,:,:),6));


%% PLOT 
% s = st;
% 
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);


nfoi=7;
cnt = 0;

for ifoi = 1:13
  
  
cnt = cnt + 1;

ss = s(:,:,:,:,ifoi);

sss(ifoi) = nanmean(ss(:));

[t1 p1]=ttest(ss(:,:,:,2),ss(:,:,:,1),'dim',3); t1(isnan(t1))=0;
[t2 p2]=ttest(ss(:,:,:,3),ss(:,:,:,1),'dim',3); t2(isnan(t1))=0;

ss_clim = squeeze(nanmean(ss,3)); ss_clim(ss_clim==inf)=nan;

clim = max([abs(min(ss_clim(:))) abs(max(ss_clim(:)))])

clim = [-clim clim];

figure_white;

subplot(1,5,1)
imagesc(nanmean(ss(:,:,:,1),3),clim); axis square; colormap(cmap)

subplot(1,5,2)
imagesc(nanmean(ss(:,:,:,2),3),clim); axis square; 

subplot(1,5,3)
imagesc(nanmean(ss(:,:,:,3),3),clim); axis square;

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_powcorr_raw_f%d_v%d.eps',ifoi,v));


figure_white;

subplot(2,2,1)
imagesc(nanmean(ss(:,:,:,2),3)-nanmean(ss(:,:,:,1),3),[-0.02 0.02]); axis square
subplot(2,2,2)
imagesc(t1,[-1 1]); axis square

subplot(2,2,3)
imagesc(nanmean(ss(:,:,:,3),3)-nanmean(ss(:,:,:,1),3),[-0.02 0.02]); axis square
subplot(2,2,4)
imagesc(t2,[-1 1]); axis square
colormap(cmap)
print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_powcorr_contrast_f%d_v%d.eps',ifoi,v));

end


%%


