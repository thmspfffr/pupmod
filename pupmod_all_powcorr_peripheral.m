%% PERMUTATION TEST OF CONNECTIVITY DIFFERENCES WITH CLEANED SIGNAL
% heart beats and blinks projected out of the FC matrices across subjects
% pupmod_all_powcorr_periphereal

% See pupmod_src_powcorr_permtest.m for statistics

clear

% -------------------------------------------------------------------------
% VERSION 01 - ignore saccades, *interpolate* blinks
% -------------------------------------------------------------------------
v = 14;
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
for icond = 1 : 2
  
  ord   = pconn_randomization;

  for m = 1:3
    for isubj = SUBJLIST
      im = find(ord(isubj,:)==m);
      
      if icond == 1
        smpdir    = dir(sprintf([sampledir_res '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        evtdir    = dir(sprintf([eventdir_res '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
      elseif icond == 2
        smpdir    = dir(sprintf([sampledir_cnt '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        evtdir    = dir(sprintf([eventdir_cnt '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
      end
      
      for iblock = 1 : length(evtdir)
        
        fprintf('Processing s%d b%d m%d ... \n',isubj,iblock,m)
        
        if icond == 1
          load([sampledir_res smpdir(iblock).name])
          load([eventdir_res evtdir(iblock).name])
          
        elseif icond == 2
          load([sampledir_cnt smpdir(iblock).name])
          load([eventdir_cnt evtdir(iblock).name])
        end
        if isubj == 26 && im == 1
          samples(:,2:4) = dat(:,1:3);
        end

        if size(samples,1) > 10000
          
          x       = abs(diff(zscore(samples(:,4))));
          [~,idx] = findpeaks(double(x>0.20),'MinPeakDistance',200);
          
          allblinks(isubj,m,icond,iblock) = length(idx);
          b(isubj,m,icond,iblock)         = length(blinks);
          
        else
          warning(sprintf('not sufficient data... s%dm%db%d',isubj,m,iblock))
          
          allblinks(isubj,m,icond,iblock) = nan;
          b(isubj,m,icond,iblock)         = nan;
          
        end
        
      end
    end
  end
end

allblinks = nanmean(allblinks(SUBJLIST,:,:,:),4);
b         = nanmean(b(SUBJLIST,:,:,:),4);

%%
for isubj = SUBJLIST
  
  for m = 1 : 3
    for ibl = 1 : 2

    try
      im = find(ord(isubj,:)==m);
      load(['~/pconn/proc/dfa/' sprintf('pconn_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v_hrv)]);
      
      hb_all(isubj,ibl,m) = par.hb;
      
       load(['~/pconn/proc/dfa/' sprintf('pconn_cnt_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v_hrv)]);

      hb_cnt(isubj,ibl,m) = par.hb;
      
    catch me
    end
    
    end
  end
end

hb_all  = hb_all(SUBJLIST,:,:);
hb_all_cnt  = hb_cnt(SUBJLIST,:,:);

hb_all(hb_all==0) = NaN;
hb_all_cnt(hb_all_cnt==0) = NaN;

hb(:,:,1)  = squeeze(nanmean(hb_all,2));
hb(:,:,2)  = squeeze(nanmean(hb_all_cnt,2));
%%
clear p1 p2

outdir = '~/pupmod/proc/conn/';
% s_fc = single(zeros(400,400,34,3,2,13));
s_fc = single(zeros(90,90,34,3,2,13));

for ifoi = 1:13
  ifoi
  
  for isubj = SUBJLIST
%     disp(isubj)
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
        
        p1(:,:,iblock) = single(powcorr);
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%        size(powcorr)
        p2(:,:,iblock) = single(powcorr);
        
        
      end
      
      s_fc(:,:,isubj,m,1,ifoi) = nanmean(p1,3);
      s_fc(:,:,isubj,m,2,ifoi) = nanmean(p2,3);
      
      clear p1 p2
      
      
    end
  end
end

s_fc = s_fc(:,:,SUBJLIST,:,:,:);

error('!')

%%

if ~exist(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v))
for icont = 1 : 2
  for im = 1 : 3
    for ifoi = 1:13
      for i = 1 : size(s_fc,1)
        fprintf('Cond %d Session %d freq %d node %d...\n',icont,im,ifoi,i)
        for j = 1 : size(s_fc,1)
      
          dat = atanh(squeeze(s_fc(i,j,:,im,icont,ifoi)))-mean(atanh(squeeze(s_fc(i,j,:,im,icont,ifoi))));
          x   = zscore(hb(:,im,icont));
          ref = x./norm(x);
         

          cleandattmp = (dat - (dat'*ref)*ref)+mean(atanh(squeeze(s_fc(i,j,:,im,icont,ifoi))));
          corr(cleandattmp,x);
          clear dat x ref
          
          dat = cleandattmp-mean(cleandattmp);
          x   = zscore(b(:,im,icont));
          ref = x./norm(x);

          cleandat(i,j,:,im,icont,ifoi) = tanh((dat - (dat'*ref)*ref)+mean(cleandattmp));

        end
      end
    end
  end
end

save(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v),'cleandat','-v7.3');
else
  load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
end

clear s_fc;

%% CONTINUE WITH PUPMOD_SRC_POWCORR_PERMTEST.M FOR STATISTICS


