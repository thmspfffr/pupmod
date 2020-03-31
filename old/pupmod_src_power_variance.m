%% pupmod_src_power_variance

clear all

% --------------------------------------------------------
% VERSION 1 - alpha0 = 0.05 (fieldtrip default)
% --------------------------------------------------------
% v         = 1;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 0.05;
% --------------------------------------------------------
% VERSION 2 - alpha0 = 0.15
% --------------------------------------------------------
% v         = 2;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 0.15;
% --------------------------------------------------------
% VERSION 3 - alpha0 = 0.3 (this is the version picked for later)
% --------------------------------------------------------
v         = 3;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid      = 'cortex_lowres';
foi_range = 2.^(2:.25:6);
REG       = 0.3;
% --------------------------------------------------------
% VERSION 4 - alpha0 = 1.0
% --------------------------------------------------------
% v         = 4;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'cortex_lowres';
% foi_range = 2.^(2:.25:6);
% REG       = 1;
% --------------------------------------------------------
% VERSION 33 - AAL and alpha0 = 0.3
% --------------------------------------------------------
% v         = 33;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% grid      = 'aal_6mm';
% foi_range = 2.^(2.75:.25:4.25); % freq range based on sign. effects 
% REG       = 0.3;
% --------------------------------------------------------

if strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
end

outdir   = '/home/tpfeffer/pupmod/proc/conn/';

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    
    if ~exist(sprintf([outdir 'pupmod_src_power_variance_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_src_power_variance_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    %
    
    fprintf('Processing s%d m%d f%d ...\n', isubj,m)
    
    for iblock = 1:2
      
      clear dat sa
      % ------------
      % Load sensor level data
      % ------------
      try
        load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('~/pupmod/proc/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))             
          if strcmp(grid,'cortex_lowres')
            powcorr(:,:,iblock,1:length(foi_range)) = nan(400,400,2,length(foi_range));
          elseif strcmp(grid,'aal_6mm')
            powcorr(:,:,iblock,1:length(foi_range)) = nan(90,90,2,length(foi_range));
          end
          continue
        else
          error('Data corrupt?')
        end
      end
      
      % ------------
      % Load source model
      % ------------
      sa        = load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      % ------------

      for ifoi = 1:length(foi_range)
        
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)
        % ------------
        % Compute cross spectrum (using wavelets)
        % ------------
        para          = [];
        para.freq     = foi_range(ifoi);
        para.fsample  = 400;           
        cs            = tp_compute_csd_wavelets(dat,para);
        
        % ------------
        % Compute spatial filter (LCMV)
        % ------------
        para          = [];
        para.reg      = REG;
        filt          = tp_beamformer(real(cs),eval(sprintf('sa.sa.L_%s',grid)),para);   
        % ------------
        % Compute variance 
        % ------------        
        wavelet = tp_mkwavelet(foi_range(ifoi),0.5,data.fsample);
        
        n_win = size(wavelet,1);
        n_shift = round(0.5*n_win);
        nseg = floor((size(dat,2)-n_win)/n_shift+1);
        
        clear datasf1
        
        for j=1:nseg
          dloc2=data.trial(:,(j-1)*n_shift+1:(j-1)*n_shift+n_win)';
          if any(isnan(dloc2(:,1)))
            warning('NaN detected')
            continue
          end
          dataf=dloc2'*KERNEL;
          datasf1(:,j)=dataf'*sa.sa.filt;
        end
        
        outp.var(:,ifoi,iblock) = var(abs(datasf1).^2,[],2);
        outp.pow(:,ifoi,iblock) = mean(abs(datasf1).^2,2);

        % ------------
        clear cs para filt
        
        
      end
    
    save(sprintf([outdir 'pupmod_src_power_variance_s%d_m%d_v%d.mat'],isubj,m,v),'outp');
    end
  end
end


error('!')

%% CLEAN NON PROCESSED FILES
outdir   = '/home/tpfeffer/pupmod/proc/conn/';

cnt = 0;
v = 6;
cnt_exist = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        ifoi
        if exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v)) && exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt_exist = cnt_exist + 1;
          
          continue
        elseif exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        elseif exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
          delete(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          cnt = cnt + 1;
        elseif ~exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && exist(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v))
          system(['touch ' outdir sprintf('pupmod_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          warning('Nothing exists')
          cnt = cnt+1;
        end
      end
    end
  end
end
cnt
%%
ord   = pconn_randomization;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:13
      for iblock = 1 : 2
        
        im = find(ord(isubj,:)==m); isubj
        
        load(sprintf([outdir 'pupmod_src_variance_s%d_m%d_v%d.mat'],isubj,m,v));
        
        var_all(:,isubj,m,ifoi,iblock) = diag(var);
        
      end
    end
  end
end

var_all = nanmean(var_all(:,SUBJLIST,:,:,:),5);

%%

for ifoi = 1 : 13
  
  h=ttest(var_all(:,:,2,ifoi),var_all(:,:,1,ifoi),'dim',2);
  n_atx(ifoi) = sum(h)./ length(h);
  
  h=ttest(var_all(:,:,3,ifoi),var_all(:,:,1,ifoi),'dim',2);
  n_dpz(ifoi) = sum(h)./ length(h);
  
end


figure; hold on

plot(n_atx,'linewidth',2,'color',[1 0.5 0.2])
plot(n_dpz,'linewidth',2,'color',[0.2 0.5 1])





%% RANGE NORMALIZE

for isubj = 1 : 28
  for ifoi = 1:21
  
  tmp = squeeze(cleandat(:,:,isubj,:,:,ifoi,:));
  min_t = min(tmp(:));
  max_t = max(tmp(:));
  
  tmp=(tmp-min_t)/(max_t-min_t);
  
  cleandat(:,:,isubj,:,:,ifoi,:)=tmp;
  
  end
end

%%

  para.alpha = 0.05;
  para.nfreq = 1:21
  
  emp = compute_altered_correlations(cleandat,para);

%%
figure_w;
subplot(3,2,1)
plot(emp.n_p_atx); axis([0 21 0 0.4])
subplot(3,2,2)
plot(emp.n_n_atx); axis([0 21 0 0.4])
subplot(3,2,3)
plot(emp.n_p_dpz); axis([0 21 0 0.4])
subplot(3,2,4)
plot(emp.n_n_dpz); axis([0 21 0 0.4])


