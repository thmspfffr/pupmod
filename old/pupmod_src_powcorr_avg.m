%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 1;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

ord = pconn_randomization;

    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));


    cleandat = cleandat./mean(cleandat,6);
    
for ifoi = 1 : 13
  
  avg_corr_rest(ifoi) = nanmean(nanmean(nanmean(cleandat(:,:,:,1,1,ifoi),3))).*1.577;
  sem_corr_rest(ifoi) = std(squeeze((nanmean(nanmean(cleandat(:,:,:,1,1,ifoi))).*1.577)))./sqrt(28);
  
  avg_corr_task(ifoi) = nanmean(nanmean(nanmean(cleandat(:,:,:,1,2,ifoi),3))).*1.577;
  sem_corr_task(ifoi) = std(squeeze((nanmean(nanmean(cleandat(:,:,:,1,2,ifoi))).*1.577)))./sqrt(28);
  
  avg_corr_rest_fish(ifoi) = nanmean(nanmean(nanmean(atanh(cleandat(:,:,:,1,1,ifoi)),3)));
  sem_corr_rest_fish(ifoi) = std(squeeze((nanmean(nanmean(atanh(cleandat(:,:,:,1,1,ifoi)))))))./sqrt(28);
  
  avg_corr_task_fish(ifoi) = nanmean(nanmean(nanmean(atanh(cleandat(:,:,:,1,2,ifoi)),3)));
  sem_corr_task_fish(ifoi) = std(squeeze((nanmean(nanmean(atanh(cleandat(:,:,:,1,2,ifoi)))))))./sqrt(28);
  
end
  


%% PLOT NUMBER OF ALTERED CONNECTIONS (IRRESPECTIVE OF SIGN)

figure;

subplot(1,2,1); hold on

shadedErrorBar(1:13,avg_corr_rest,sem_corr_rest)
shadedErrorBar(1:13,avg_corr_task,sem_corr_task)

axis([0 14 0 4])
set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

xlabel('Frequency [Hz]'); ylabel('Number of alt. corr. [%]')
axis square
box off
% 
% subplot(1,2,2); hold on
% 
% shadedErrorBar(1:13,avg_corr_rest_fish,sem_corr_rest_fish)
% shadedErrorBar(1:13,avg_corr_task_fish,sem_corr_task_fish)
% 
% axis([0 14 0 4])
% set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
% 
% xlabel('Frequency [Hz]'); ylabel('Number of alt. corr. [%]')
% axis square
% box off

print(gcf,'-depsc2',sprintf('~/pupmod/plots/pupmod_src_powcorr_avgconn.pdf'))

%% PLOT
