%% OBTAIN ALL PERIPHERAL SIGNALS
% Regress out in next script: pupmod_all_regressartifacts.m
% See pupmod_src_powcorr_permtest.m for statistics

clear

% -------------------------------------------------------------------------
% VERSION 01 - ignore saccades, *interpolate* blinks
% -------------------------------------------------------------------------
v = 1;
v_pup = 2;
v_hrv = 1;
% -------------------------------------------------------------------------

tp_addpaths

sampledir_cnt   = '/home/tpfeffer/pconn_cnt/proc/';
eventdir_cnt   = '/home/tpfeffer/pconn_cnt/proc/';
outdir      = '/home/tpfeffer/pconn_cnt/proc/';

sampledir_res   = '/home/tpfeffer/pconn/proc/pup/';
eventdir_res    = '/home/tpfeffer/pconn/proc/pup/';


SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pconn/matlab
%%
% missing blinks are taken from pconn_find_blinks.m

for icond = 1 : 2
  
  ord   = pconn_randomization;
  
  for m = 1:3
    for isubj = SUBJLIST
      im = find(ord(isubj,:)==m);
      
      clear smpdir evtdir
      
      if icond == 1
        smpdir    = dir(sprintf([sampledir_res '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        evtdir    = dir(sprintf([eventdir_res '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
      elseif icond == 2
        smpdir    = dir(sprintf([sampledir_cnt '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        evtdir    = dir(sprintf([eventdir_cnt '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
      end
      
      if length(evtdir) == 0
        if icond == 1
          load(['~/pconn/proc/preproc/' sprintf('pconn_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,1,2)])
          b(isubj,m,icond,1) = blinks;
          allblinks(isubj,m,icond,1) =blinks; clear blinks
          load(['~/pconn/proc/preproc/' sprintf('pconn_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,2,2)])
          b(isubj,m,icond,2) = blinks;
          allblinks(isubj,m,icond,2) =blinks; clear blinks
        else
          load(['~/pconn/proc/preproc/' sprintf('pconn_cnt_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,1,1)])
          b(isubj,m,icond,1) = blinks;
          allblinks(isubj,m,icond,1) =blinks; clear blinks
          load(['~/pconn/proc/preproc/' sprintf('pconn_cnt_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,2,1)])
          b(isubj,m,icond,2) = blinks;
          allblinks(isubj,m,icond,2) =blinks; clear blinks
        end
        warning(sprintf('No files: s%d m%d c%d', isubj,im,icond))
        continue
      end
      
      for iblock = 1 : length(evtdir)
        
        %         clear evtdir
        fprintf('Processing s%d b%d m%d ... \n',isubj,iblock,m)
        
        load([smpdir(iblock).folder '/' smpdir(iblock).name])
        load([evtdir(iblock).folder '/' evtdir(iblock).name])
        blinks1 = blinks;
        
        if icond == 1
          if isubj < 10
            block = str2num(evtdir(iblock).name(22));
          else
            block = str2num(evtdir(iblock).name(23));
          end
        else
          if isubj < 7
            block = str2num(evtdir(iblock).name(26));
          else
            block = str2num(evtdir(iblock).name(27));
          end
        end
        
        if length(evtdir) == 1
          
          if block == 1
            if icond == 1
              load(['~/pconn/proc/preproc/' sprintf('pconn_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,2,2)])
            else
              load(['~/pconn/proc/preproc/' sprintf('pconn_cnt_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,2,1)])
            end
            
            b(isubj,m,icond,2) = blinks;
            allblinks(isubj,m,icond,2) = blinks; clear blinks
            
            warning(sprintf('Missing block: s%d m%d b%d c%d', isubj,im,2,icond))
          else
            if icond == 1
              load(['~/pconn/proc/preproc/' sprintf('pconn_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,1,2)])
            else
              load(['~/pconn/proc/preproc/' sprintf('pconn_cnt_find_blinks_s%d_m%d_b%d_v%d.mat',isubj,im,1,1)])
            end
            b(isubj,m,icond,1) = blinks;
            allblinks(isubj,m,icond,1) = blinks; clear blinks
            warning(sprintf('Missing block: s%d m%d b%d c%d', isubj,im,1,icond))
          end
        end
        
        if isubj == 26 && im == 1
          samples(:,2:4) = dat(:,1:3);
        end
        
        if size(samples,1) > 10000
          samples = samples(10000:end,:);
          x = highpass(samples(:,4),50,1000);
          x(x<50) = 0;
          %           x       = -(zscore(diff(samples(:,4))));
          [~,idx] = findpeaks(double(x),'MinPeakDistance',500);
          
          allblinks(isubj,m,icond,block) = length(idx);
          
        else
          
          allblinks(isubj,m,icond,block) = nan;
        end
        
        try
          b(isubj,m,icond,block)         = sum(blinks1(:,2)>samples(1,1)); clear blinks1
          
        catch me
          b(isubj,m,icond,block)         = nan;
        end
        
      end
    end
  end
end

allblinks = allblinks(SUBJLIST,:,:,1:2);
b         = b(SUBJLIST,:,:,1:2);

%% MISSING HEART RATE INFO TAKEN FROM ICA COMPONENTS
% see pconn_hrv_fromica.m for details (same adaptivethresholding procedure
% is applied)

clear hb_all hb_cnt

for isubj = SUBJLIST
  
  for m = 1 : 3
    for ibl = 1 : 2
      
      try
        im = find(ord(isubj,:)==m);
        load(['~/pconn/proc/dfa/' sprintf('pconn_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v_hrv)]);
        hb_all(isubj,m,ibl) = par.hb;
        
      catch me
        try
          load(sprintf('~/pconn/proc/pconn_hrv_fromica_s%d_m%d_b%d.mat',isubj,im,ibl),'par')
          hb_all(isubj,m,ibl) = par;
        catch me
          hb_all(isubj,m,ibl) = nan;
          warning(sprintf('Did not find REST data: s%d m%d ibl%d',isubj,im,ibl))
          
        end
      end
      
      try
        
        load(['~/pconn/proc/dfa/' sprintf('pconn_cnt_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v_hrv)]);
        
        hb_cnt(isubj,m,ibl) = par.hb;
        
      catch me
        try
          load(sprintf('~/pconn_cnt/proc/pconn_cnt_hrv_fromica_s%d_m%d_b%d.mat',isubj,im,ibl),'par')
          hb_cnt(isubj,m,ibl) = par;
        catch me
          hb_cnt(isubj,m,ibl) = nan;
          warning(sprintf('Did not find TASK data: s%d m%d ibl%d',isubj,im,ibl))
        end
      end
      
      
    end
  end
end

clear hb

hb_all  = hb_all(SUBJLIST,:,:);
hb_all_cnt  = hb_cnt(SUBJLIST,:,:);

hb_all(hb_all==0) = NaN;
hb_all_cnt(hb_all_cnt==0) = NaN;

hb(:,:,1,:)  = squeeze(hb_all);
hb(:,:,2,:)  = squeeze(hb_all_cnt);

%% OTHER ARTIFACTS
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for ibl = 1 : 2
      try
        load(['/home/tpfeffer/pconn/proc/preproc/' sprintf('pconn_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,1)])
        artcnt(isubj,m,1,ibl) = size(art,1);
        artlen(isubj,m,1,ibl) = mean(art(:,2)-art(:,1));
      catch me
        try
          load(['/home/tpfeffer/pconn/proc/preproc/' sprintf('pconn_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,2)])
          artcnt(isubj,m,1,ibl) = size(art,1);
          artlen(isubj,m,1,ibl) = mean(art(:,2)-art(:,1));
        catch
          artcnt(isubj,m,1,ibl) = nan;
          artlen(isubj,m,1,ibl) = nan;
        end
      end
      
    end
  end
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for ibl = 1 : 2
      try
        load(['/home/tpfeffer/pconn_cnt/proc/preproc/' sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,1)])
        artcnt(isubj,m,2,ibl) = size(art,1);
        artlen(isubj,m,2,ibl) = mean(art(:,2)-art(:,1));
      catch me
        try
          load(['/home/tpfeffer/pconn_cnt/proc/preproc/' sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,1)])
          artcnt(isubj,m,2,ibl) = size(art,1);
          artlen(isubj,m,2,ibl) = mean(art(:,2)-art(:,1));
        catch
          artcnt(isubj,m,2,ibl) = nan;
          artlen(isubj,m,2,ibl) = nan;
        end
      end
    end
  end
end

artcnt= artcnt(SUBJLIST,:,:,:);

%%
par = [];
par.heartrate = hb;
par.blinks = b;
par.muscles = artcnt;

save([outdir sprintf('pupmod_all_powcorr_peripheral.mat')], 'par')

