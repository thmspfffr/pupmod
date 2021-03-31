%% pupmopd_src_peakfreq
% compute powerspectra in source space using LCMV (broadband)
% for later application of fooof (in python)
% 10-02-2021

clear
restoredefaultpath

% -------------------------
% VERSION 33 - AAL and alpha0 = 0.3
% --------------------------------------------------------
v         = 33;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid      = 'aal_6mm';
foi_range = 2.^(2:(1/4):6); % freq range based on sign. effects
REG       = 0.3;
% --------------------------------------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pupmod/matlab/

ft_defaults
addpath ~/pconn/matlab/
outdir = '~/pupmod/proc/src/';

ord    = pconn_randomization;

if strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
end

%%
% -------------------------
for isubj = SUBJLIST
  
  for m = 1 : 3
    
    for iblock = 1:2
      %
      fn = sprintf('pupmod_rest_peakfreq_s%d_m%d_b%d_v%d',isubj,m,iblock,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
      %
      fprintf('Processing subj%d block%d ...\n',isubj,iblock);
      
      try
        % load cleaned meg data
        load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
      catch me
        continue
      end
      
      load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      lf = sa.L_aal_6mm;
      
      csd = [];
      for ifreq=1:length(foi_range)
        ifreq
        
        % -------------------------------
        % compute CSD
        % -------------------------------
        para          = [];
        para.freq     = foi_range(ifreq);
        para.fsample  = 400;
        para.overlap  = 0.5;
        csd(:,:,ifreq) = tp_compute_csd_wavelets(dat,para);
        
      end
      
      % get broad band CSD
      csd = nanmean(csd,3);
      
      dat = dat(:,~isnan(dat(1,:)));
        % -------------------------------
        % beamforming
        % -------------------------------
        para          = [];
        para.reg      = 0.05;
        [filt,outp.pow(:,ifreq)] = tp_beamformer(real(csd),lf,para);
        % --------------
        
        for ireg = 1 : size(sa.aal_centr,1)
          idx=find(sa.aal_label==ireg);       
          src_dat = dat'*filt(:,idx);
          [pxx,outp.fxx]=pwelch(src_dat,hanning(8000),4000,3:0.05:40,400,'power');
          
          outp.pxx(:,ireg) = nanmean(pxx,2);
        end

      
      
      save([outdir fn '.mat'],'outp')
      tp_parallel(fn,outdir,0)
      
      clear outp
      end
    end
  end



% DO SAME FOR TASK DATA
for isubj = SUBJLIST
  
  for m = 1 : 3
    
    for iblock = 1:2
      %
      fn = sprintf('pupmod_task_peakfreq_s%d_m%d_b%d_v%d',isubj,m,iblock,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
      %
      fprintf('Processing subj%d block%d ...\n',isubj,iblock);
      
      try
        % load cleaned meg data
        load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
      catch me
        continue
      end
      
      load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid));
      lf = sa.L_aal_6mm;
      
      csd = [];
      for ifreq=1:length(foi_range)
        ifreq
        
        % -------------------------------
        % compute CSD
        % -------------------------------
        para          = [];
        para.freq     = foi_range(ifreq);
        para.fsample  = 400;
        para.overlap  = 0.5;
        csd(:,:,ifreq) = tp_compute_csd_wavelets(dat,para);
        
      end
      
      % get broad band CSD
      csd = nanmean(csd,3);
      dat = dat(:,~isnan(dat(1,:)));
        % -------------------------------
        % beamforming
        % -------------------------------
        para          = [];
        para.reg      = 0.05;
        [filt,outp.pow(:,ifreq)] = tp_beamformer(real(csd),lf,para);
        % --------------
        clear outp
        for ireg = 1 : size(sa.aal_centr,1)
          idx=find(sa.aal_label==ireg);       
          src_dat = dat'*filt(:,idx);
          [pxx,outp.fxx]=pwelch(src_dat,hanning(8000),4000,3:0.05:40,400,'power');
          
          outp.pxx(:,ireg) = nanmean(pxx,2);
        end

      
      
      save([outdir fn '.mat'],'outp')
      tp_parallel(fn,outdir,0)
      
      
      
    end
  end
end


exit

%%
fc_tmp = pupmod_loadpowcorr(33,SUBJLIST,1);
fc_tmp=nanmean(fc_tmp(:,:,:,:,:,3:6),6);

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

for isubj = 1 : 28
  for im = 1 : 3
    for icond = 1 : 2
      tmp = squeeze(fc_tmp(:,:,isubj,im,icond));
      tmp = tp_match_aal(para,tmp); 
      fc(:,:,isubj,im,icond) = tmp(include_bcn,include_bcn);
    end
  end
end
 
[h_atx,p]=ttest(squeeze(nanmean(fc(:,:,:,2,2),1)),squeeze(nanmean(fc(:,:,:,1,2),1)),'dim',2)
[h_dpz,p]=ttest(squeeze(nanmean(fc(:,:,:,3,1),1)),squeeze(nanmean(fc(:,:,:,1,1),1)),'dim',2)



      
%% LOAD FOOOF RESULST (from pupmod_src_fooof.py)
para          = [];
para.transfer = 'to_bcn';
para.N        = 90;
[~,trans,lab]=tp_match_aal(para,rand(90,90));
load(sprintf('~/pupmod/proc/src/pupmod_src_pow_taskmodulationindex_v%d.mat',33))

  
k = 1 : 90;

% exclude subcortical regions
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));


v = 33;

fxx = 3:0.05:40;
alpha_idx = fxx>6 & fxx<14;


SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
ord    = pconn_randomization;
 
clear all_pow all_pow2

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      
      im = find(ord(isubj,:)==m);
      
      try 
        load(sprintf('~/pupmod/proc/src/pupmod_rest_fooof_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
        all_pow(:,:,isubj,m,1,iblock) = only_gauss(:,trans(:,2));
        all_pow2(:,:,isubj,m,1,iblock) = full_gauss(:,trans(:,2));
        all_aper(:,:,isubj,m,1,iblock) = aper(:,:);
      catch me
        all_pow(:,:,isubj,m,1,iblock) = nan(741,90);
        all_pow2(:,:,isubj,m,1,iblock) = nan(741,90);
        all_aper(:,:,isubj,m,1,iblock) = nan(2,91);
      end
      
      try
        load(sprintf('~/pupmod/proc/src/pupmod_task_fooof_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
        all_pow(:,:,isubj,m,2,iblock) = only_gauss(:,trans(:,2));
        all_pow2(:,:,isubj,m,2,iblock) = full_gauss(:,trans(:,2));
        all_aper(:,:,isubj,m,2,iblock) = aper(:,:);
      catch me
        all_pow(:,:,isubj,m,2,iblock) = nan(741,90);
        all_pow2(:,:,isubj,m,2,iblock) = nan(741,90);
        all_aper(:,:,isubj,m,2,iblock) = nan(2,91);
      end
      
      
      
    end
  end
end

all_pow = squeeze(all_pow(:,include_bcn,SUBJLIST,:,:,:));
all_pow2 = squeeze(all_pow2(:,include_bcn,SUBJLIST,:,:,:));
% all_aper = squeeze(all_aper(:,include_bcn,SUBJLIST,:,:,:));

all_pow = all_pow(:,:,:,:,:,:);
alpha_idx=find(alpha_idx);

%% LOAD EMPIRICAL POWER
% and remove 1/f component (aper_all)

clear emp_pow
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      
      im = find(ord(isubj,:)==m);
      
      
      try
        load(sprintf('~/pupmod/proc/src/pupmod_rest_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
        aper_comp = repmat(all_aper(1,:,isubj,m,1,iblock),[741 1])'-all_aper(2,:,isubj,m,1,iblock)'.*log10((3:0.05:40));
        emp_pow(:,:,isubj,m,1,iblock) = log10(outp.pxx(:,trans(:,2)))-aper_comp(trans(:,2),:)';
        

      catch me
        emp_pow(:,:,isubj,m,1,iblock) = nan(741,90);
      end
      
      try
        load(sprintf('~/pupmod/proc/src/pupmod_task_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
        aper_comp = repmat(all_aper(1,:,isubj,m,2,iblock),[741 1])'-all_aper(2,:,isubj,m,2,iblock)'.*log10((3:0.05:40));
        emp_pow(:,:,isubj,m,2,iblock) = log10(outp.pxx(:,trans(:,2)))-aper_comp(trans(:,2),:)';
      catch me
        emp_pow(:,:,isubj,m,2,iblock) = nan(741,90);
      end
      
    end
  end
end

emp_pow = emp_pow(:,include_bcn,SUBJLIST,:,:,:);

%% FIND PEAK

clear pf

for isubj = 1 : 28
  isubj
  for icond = 1 : 2
    for im =1 :3
      for ireg = 1 : size(all_pow,2)
        for iblock = 1 : 2
          
          [mag,loc]=findpeaks(squeeze(nanmean(all_pow(alpha_idx,ireg,isubj,im,icond,iblock),6)));
          if ~isempty(mag)
%             pf(ireg,isubj,im,icond,iblock) = fxx(alpha_idx(loc(sum((mag./mag')>=1)==length(mag))));
            pf(ireg,isubj,im,icond,iblock) = sum((mag./sum(mag)) .* fxx(alpha_idx(loc))',1);
          else
            pf(ireg,isubj,im,icond,iblock) = nan;
          end
          
        end
      end
    end
  end
end

% idx=find(sum(~isnan(pf(:,:,1,1)),2)==28)

pf = squeeze(nanmean(pf,1));

figure_w;
subplot(2,3,1); 
scatter(pf(:,1,1,1),pf(:,1,1,2),40,'markeredgecolor','w','markerfacecolor','k')
tp_editplots; axis equal square; xlabel('Peak freq. (Block 1) [Hz]'); ylabel('Peak freq. (Block 2) [Hz]')
axis([8 12 8 12]); lsline
[r,p]=corr(pf(:,1,1,1),pf(:,1,1,2));
text(8.3,11.6,sprintf('r = %.3f | p = %.3f',r,p),'fontsize',6)

subplot(2,3,2); 
scatter(pf(:,2,1,1),pf(:,2,1,2),40,'markeredgecolor','w','markerfacecolor','k')
tp_editplots; axis equal square; xlabel('Peak freq. (Block 1) [Hz]'); ylabel('Peak freq. (Block 2) [Hz]')
axis([8 12 8 12]); lsline
[r,p]=corr(pf(:,2,1,1),pf(:,2,1,2));
text(8.3,11.6,sprintf('r = %.3f | p = %.3f',r,p),'fontsize',6)

subplot(2,3,3); 
scatter(pf(:,3,1,1),pf(:,3,1,2),40,'markeredgecolor','w','markerfacecolor','k')
tp_editplots; axis equal square; xlabel('Peak freq. (Block 1) [Hz]'); ylabel('Peak freq. (Block 2) [Hz]')
axis([8 12 8 12]); lsline
[r,p]=corr(pf(:,3,1,1),pf(:,3,1,2));
text(8.3,11.6,sprintf('r = %.3f | p = %.3f',r,p),'fontsize',6)

pf=nanmean(pf,4);

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_peakfreq_controls_v%d.pdf',v))


%%
marker = 4;
randnumb = (rand(28,1)-0.5)/3;

greys = cbrewer('seq', 'Greys', 120,'pchip'); greys=greys(32:2:end-32,:)

figure_w

subplot(3,2,2); hold on

plot(fxx,nanmean(nanmean(nanmean(all_pow(:,:,:,1,1,:),2),3),6),'k:')
plot(fxx,nanmean(nanmean(nanmean(all_pow(:,:,:,1,2,:),2),3),6),'k')

plot(fxx,nanmean(nanmean(nanmean(emp_pow(:,:,:,1,1,:),2),3),6),'color',[0.5 0.5 0.5],'linestyle',':')
plot(fxx,nanmean(nanmean(nanmean(emp_pow(:,:,:,1,2,:),2),3),6),'color',[0.5 0.5 0.5],'linestyle','-')

[~,k1]=max(nanmean(nanmean(nanmean(all_pow(alpha_idx,:,:,1,1,:),2),3),6));
[~,k2]=max(nanmean(nanmean(nanmean(all_pow(alpha_idx,:,:,1,2,:),2),3),6));

line([fxx(alpha_idx(k1)) fxx(alpha_idx(k1))],[0 0.5],'color','k','linestyle',':')
line([fxx(alpha_idx(k2)) fxx(alpha_idx(k2))],[0 0.5],'color','k','linestyle','-')

tp_editplots; axis square
axis([3 16 0 0.4])

subplot(3,2,1); hold on

for i = 1 : 28
  plot(ones(1,1)+randnumb(i),pf(i,1,1),'o','markersize',marker,'markeredgecolor','none','markerfacecolor',greys(i,:))
  plot(ones(1,1)*2+randnumb(i),pf(i,1,2),'o','markersize',marker,'markeredgecolor','none','markerfacecolor',greys(i,:))
  line([1+randnumb(i) 2+randnumb(i)],[pf(i,1,1) pf(i,1,2)],'color',greys(i,:))
end

line([0.6 1.4],[nanmean(pf(:,1,1)) nanmean(pf(:,1,1))],'color','k','linewidth',2)
line([1.6 2.4],[nanmean(pf(:,1,2)) nanmean(pf(:,1,2))],'color','r','linewidth',2)

axis([0 3 6 14]); axis square; tp_editplots; ylabel('Peak frequency'); xlabel('Rest / Task')

p = tp_permtest(pf(:,1,2),pf(:,1,1),100000,0)
  
print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_peakfreq_v%d.pdf',v))
%   
% 
%%
diff_pf = pf(:,1,2)-pf(:,1,1);
fc = pupmod_loadpowcorr(33,SUBJLIST,1);

tmp=squeeze(nanmean(nanmean(fc(:,:,:,2,2,3:6),2),6)-nanmean(nanmean(fc(:,:,:,1,2,3:6),2),6))

clear diff_atx
for isubj = 1 : 28
  tmp1  = tmp(trans(:,2),isubj);
  diff_atx(:,isubj)=tmp1(include_bcn,:);
end

diff_atx = nanmean(diff_atx,1);
% k = 1 : 90;
% 
% % exclude subcortical regions
% exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
% include_bcn = find(~ismember(k,exclude_bcn));
% 
% % load SC matrix, exclude subcortical regions
% load ~/sc90.mat
% SC = SC(include_bcn,include_bcn);
%   
%   
%   




