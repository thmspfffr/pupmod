% pupmod_all_pac


v               = 1;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

low_freqs = 3:1:20;

high_freqs = 40:2:80;

outdir   = '/home/tpfeffer/pupmod/proc/pac/';

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
addpath ~/pconn/matlab/

%%


for isubj = SUBJLIST
  for m = 1 : 3
    for ilowf = 1:length(low_freqs)
      
      if ~exist(sprintf([outdir 'pupmod_all_pac_s%d_m%d_lf%d_v%d_processing.txt'],isubj,m,ilowf,v))
        system(['touch ' outdir sprintf('pupmod_all_pac_s%d_m%d_lf%d_v%d_processing.txt',isubj,m,ilowf,v)]);
      else
        continue
      end
      
      
      for iblock = 1:2
        
        fprintf('Loading MEG data ...\n');
        
        try
          load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          if ~exist(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
            
            mi = nan(275,length(high_freqs));
            
            save(sprintf('~/pupmod/proc/pac/pupmod_all_pac_s%d_m%d_lf%d_v%d.mat',isubj,m,ilowf,v),'mi')
            continue
          else
            error('Data corrupt?')
          end
        end
        
        if isubj >= 32
          load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
        elseif isubj < 4
          load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
        elseif isubj == 17
          load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
        else
          load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
        end
        
        for ihighf = 1:length(high_freqs)
          
          hf = [high_freqs(ihighf)-low_freqs(ilowf) high_freqs(ihighf)+low_freqs(ilowf)];
          lf = [low_freqs(ilowf)-1 low_freqs(ilowf)+1];
          
          fprintf('Processing s%d m%d lf%d hf%d ...\n', isubj,m,ilowf,ihighf)
          
          dat_low = angle(hilbert(ft_preproc_bandpassfilter(data.trial{1},400,lf)));
          dat_hi  = abs(hilbert(ft_preproc_bandpassfilter(data.trial{1},400,hf)));
          
          fprintf('Computing PAC...\n')
          
          phase_degrees = dat_low*360./2./pi; % Phases in degrees
          clear dat_low
          n_bins = 18;
          n_chan = size(dat_hi,1);
          
          % Bining the phases
          step_length = 360/n_bins;
          phase_bins = -180:step_length:180;
          [~,phase_bins_ind] = histc(phase_degrees,phase_bins);
          clear phase_degrees
          
          % Averaging amplitude time series over phase bins
          amplitude_bins = zeros(n_chan,n_bins);
          n_chan = size(dat_hi,1);
          
          for ichan = 1:n_chan
            for bin = 1:n_bins
              amplitude_bins(ichan,bin) = mean(dat_hi(ichan,phase_bins_ind(ichan,:)==bin),2);
            end
          end
          
          clear dat_hi
          % Normalize amplitudes
          P = amplitude_bins./repmat(sum(amplitude_bins,2),1,n_bins);
          % Compute modulation index
          tmp_mi(:,ihighf,iblock) = 1+sum(P.*log(P),2)./log(n_bins);
          
        end
      end
      
      imi = mean(tmp_mi,3); clear tmp_mi
      
      for i = 1 : size(imi,2)
         mi(:,i) = pconn_sens_interp274(idx,imi(:,i));
      end
      
      clear imi
      
      save(sprintf('~/pupmod/proc/pac/pupmod_all_pac_s%d_m%d_lf%d_v%d.mat',isubj,m,ilowf,v),'mi')
      
      clear mi
      
    end
  end
end

%% LOAD DATA
addpath ~/pconn/matlab/
ord   = pconn_randomization;

for isubj =4:34
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for ilowf = 1:length(low_freqs)
      
     
      try
        load(sprintf('~/pupmod/proc/pac/pupmod_all_pac_s%d_m%d_lf%d_v%d.mat',isubj,m,ilowf,v))
        all_mi(:,:,isubj,m,ilowf) = mi;        
      catch me
        warning('...')
        all_mi(:,:,isubj,m,ilowf) = nan(274,21);
      end
            
    end
  end
end

all_mi   = all_mi(:,:,SUBJLIST,:,:);
%%
% all_mi(:,:,27,:,:) = [];
load redblue.mat
cmap = redblue;

subplot(1,1,1); 
[~,~,~,t]=ttest(squeeze(nanmean(all_mi(:,:,:,3,:),1)),squeeze(nanmean(all_mi(:,:,:,1,:),1)),'dim',2);
imagesc(squeeze(t.tstat),[0 2.5])
set(gca,'tickdir','out','xtick',1:length(low_freqs),'xticklabel',num2cell([low_freqs]))
set(gca,'tickdir','out','ytick',1:length(high_freqs),'yticklabel',num2cell([high_freqs]))
set(gca,'ydir','normal')
colormap(redblue)