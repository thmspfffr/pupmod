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

outdir   = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/


%%
clear powcorr_all
v = 10;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/
  
ord = pconn_randomization;

for ifoi = 1:13
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
       
        powcorr_all(:,:,isubj,m,1,ifoi,iblock) =  powcorr; clear powcorr
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));

        powcorr_all(:,:,isubj,m,2,ifoi,iblock) =  powcorr; clear powcorr
              
      end
    end
  end
end

powcorr_all = squeeze(nanmean(powcorr_all(:,:,SUBJLIST,:,:,:,:),7));


%% PLOT 


for ifoi = 1 : 13
  
  fcmat_rest = squeeze(nanmean(powcorr_all(:,:,:,1,2,ifoi),3));
  fcmat_task = squeeze(nanmean(powcorr_all(:,:,:,2,2,ifoi),3));

  fcmat_diff = squeeze(nanmean(powcorr_all(:,:,:,1,2,ifoi),3))-squeeze(nanmean(powcorr_all(:,:,:,1,1,ifoi),3));
  
  m_fc_diff(ifoi) = mean(fcmat_diff(ask_find(triu(ones(90,90),1))));
  m_fc_rest(ifoi) = mean(fcmat_rest(find(triu(ones(90,90),1))));
  m_fc_task(ifoi) = mean(fcmat_task(find(triu(ones(90,90),1))));
  

  
end


figure;
subplot(1,2,1); hold on
plot(m_fc_rest,'r','linewidth',3)
plot(m_fc_task,'color',[0.2 0.5 1],'linewidth',3)

axis square; xlabel('Frequency [Hz]'); ylabel('Correlation');
axis([0 14 0.015 0.035 ]); set(gca,'xtick',[1 5 9 13],'xticklabel',[2 8 32 128])



%%


