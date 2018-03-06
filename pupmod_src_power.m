%% pupmod_src_power

% Script computes center frequency for alpha/beta band in order
% to investigate a potential frequency shift due to neuromodulation
% as predicted by simulations of a WC model (email Adrian 19/02/18)

% v1, 26/02/18

clear

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v               = 1;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
%       f = 8:0.1:14;

% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v               = 2;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% f = [7 13; 14 22];
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v               = 3;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% f               = [7 13]; % relevant for beamformer csd
% segleng         = 10000;
% --------------------------------------------------------
% VERSION 4 - GAMMA
% --------------------------------------------------------
% v               = 4;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% f               = [50 100]; % relevant for beamformer csd
% segleng         = 400;
% --------------------------------------------------------
% VERSION 5 - BROADBAND
% --------------------------------------------------------
v               = 4;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'vtpm_6mm';
f               = [2 100]; % relevant for beamformer csd
segleng         = 400;
% --------------------------------------------------------


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pupmod/proc/pow/';
addpath /home/tpfeffer/pconn/matlab/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 400;

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid,'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'vtpm_4mm')
  v_grid = 10;
elseif strcmp(allpara.grid,'vtpm_6mm')
  v_grid = 11;
end


%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % %
    if ~exist(sprintf([outdir 'pupmod_src_power_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_src_power_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    %
    fprintf('Processing s%d m%d f%d ...\n', isubj,m)
    
    for iblock = 1:2
      
      fprintf('Loading MEG data ...\n');
      
      
      load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      [dat] = megdata2mydata(data); clear data
      
      pars      = [];
      pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      load(pars.sa);
      
      pars = [];
      pars.grid      = allpara.grid;
      
      fsample = 400;
      all_freq = 0:fsample/segleng:fsample/2;
      
      segshift = segleng / 2;
      maxfreqbin = fsample/2;
      epleng    = size(dat,1);
      f_select = find(all_freq>=f(1) & all_freq<=f(2));
      
      cs = data2cs_event(dat,segleng,segshift,epleng,maxfreqbin);
      
      tmp = regexp(pars.grid,'_','split');
      
      para.reg = 0.05;
      para.iscs = 1;
      
      eval(sprintf('pos = sa.grid_%s_indi;',sprintf('%s%s',tmp{1},tmp{2})));
      eval(sprintf('filt = pconn_beamformer(nanmean(cs(:,:,f_select),3),sa.L_%s,para);',pars.grid));
        
      for ifoi = 1 : maxfreqbin
        ifoi
        for iaal = 1 : 90
          cs_src    = filt(:,sa.aal_label==iaal)'*squeeze(cs(:,:,ifoi))*filt(:,sa.aal_label==iaal);
          pow(iaal,ifoi,1) = mean(real(diag(cs_src)));
        end
      end
      
      clear sa dat data cs cs_src 
      try
        
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        [dat] = megdata2mydata(data); clear data

        pars      = [];
        pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        load(pars.sa);
        
        pars = [];
        pars.grid      = allpara.grid;
        
        fsample = 400;
        all_freq = 0:fsample/segleng:fsample/2;
        
        segshift = segleng / 2;
        epleng    = size(dat,1);
        f_select = find(all_freq>=f(1) & all_freq<=f(2));
        
        cs = data2cs_event(dat,segleng,segshift,epleng,maxfreqbin);
        
        tmp = regexp(pars.grid,'_','split');
        
        para.reg = 0.05;
        para.iscs = 1;
        
        eval(sprintf('pos = sa.grid_%s_indi;',sprintf('%s%s',tmp{1},tmp{2})));
        eval(sprintf('filt = pconn_beamformer(nanmean(cs(:,:,f_select),3),sa.L_%s,para);',pars.grid));
        
        
        for ifoi = 1 : maxfreqbin
          ifoi
          for iaal = 1 : 90
            cs_src    = filt(:,sa.aal_label==iaal)'*squeeze(cs(:,:,ifoi))*filt(:,sa.aal_label==iaal);
            pow(iaal,ifoi,2) = mean(real(diag(cs_src)));
          end
        end
      catch me
        pow(:,:,2) = nan(size(pow,1),size(pow,2));
        save(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],isubj,m,iblock,v),'all_freq','pow','f_select');

      end
      
      save(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],isubj,m,iblock,v),'all_freq','pow','f_select');
      
      
    end
  end
end

error('!')

tp_addpaths
%%
clear fc_all pow_all pow_all_res pow_all_tsk  fc_all_res fc_all_tsk
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25   30 31 32  34];
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
ord = pconn_randomization;
peak_range = [7 13];
cmap = cbrewer('qual', 'Paired', 12,'pchip');
cmap = cmap([2 6 8 1 5 7],:);
set(groot,'defaultaxescolororder',cmap);

v = 3;

load(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],5,1,1,v));


for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 :2
      
      load(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,v));
      
      pow = pow(:,find(all_freq>=2,1,'first'):end,:);
      all_freq = all_freq(all_freq>=2);
      
      f_select = find(all_freq>=peak_range(1) & all_freq<=peak_range(2));
      
      para.f       = all_freq;
      para.fselect = f_select;
      para.method  = 'com';
      para.win     = 5;
      para.detrend = 0;
      
      if para.detrend == 1      
        tmp_pow = log10(squeeze(nanmean(pow)));
        pow_detr = exp(detrend(tmp_pow));     
      else
        pow_detr = squeeze(nanmean(pow));
      end
       
      fc_all_res(isubj,m,iblock) = tp_peakfreq(pow_detr(:,1),para);      
      
      try
        fc_all_tsk(isubj,m,iblock) = tp_peakfreq(pow_detr(:,2),para);
      catch me
        fc_all_tsk(isubj,m,iblock) = nan;
      end
%       [~,i]=max(nanmean(pow(:,f_select,1),1));
%       fc_all_res(isubj,m,iblock) = all_freq(f_select(i));
%       pos_res(isubj,m,iblock) = i;
%       [~,i]=max(nanmean(pow(:,f_select,2),1));
%       fc_all_tsk(isubj,m,iblock) = all_freq(f_select(i));
%       pos_tsk(isubj,m,iblock) = i;
      
      pow_all_res(:,isubj,m,iblock) = nanmean(pow(:,:,1),1);
      pow_all_tsk(:,isubj,m,iblock) = nanmean(pow(:,:,2),1);
      
    end
  end
end

fc_all_res  = fc_all_res(SUBJLIST,:,:);
mean_fc_res = nanmean(fc_all_res,3);
pow_all_res = nanmean(pow_all_res(:,SUBJLIST,:,:),4);

fc_all_tsk  = fc_all_tsk(SUBJLIST,:,:);
mean_fc_tsk = nanmean(fc_all_tsk,3);
pow_all_tsk = nanmean(pow_all_tsk(:,SUBJLIST,:,:),4);

% fc_all(1) = f*mean(pow_all(:,:,1),2) / sum(mean(pow_all(:,:,1),2));
% fc_all(2) = f*mean(pow_all(:,:,2),2) / sum(mean(pow_all(:,:,2),2));
% fc_all(3) = f*mean(pow_all(:,:,3),2) / sum(mean(pow_all(:,:,3),2));

figure; hold on
subplot(3,2,[1 3]); hold on

mp_res = squeeze(nanmean(pow_all_res,2));
mp_tsk = squeeze(nanmean(pow_all_tsk,2));

mp_res = detrend(mp_res);
mp_tsk = detrend(mp_tsk);
plot(all_freq(1:length(mp_res)),mp_res,'linewidth',2)
plot(all_freq(1:length(mp_tsk)),mp_tsk,'linewidth',2)
axis square

% m_res = squeeze(nanmean(mean_fc_res));
% m_tsk = squeeze(nanmean(mean_fc_tsk));
      para.method  = 'peak';

for i = 1 : 3
  m_res(i) = tp_peakfreq(squeeze(nanmean(pow_all_res(:,:,i),2)),para);
  m_tsk(i) = tp_peakfreq(squeeze(nanmean(pow_all_tsk(:,:,i),2)),para);
end

% [~,i_res]=max(mean(mean(pow_all_res(f_select,:,:,:),4),2))
% [~,i_tsk]=max(mean(mean(pow_all_tsk(f_select,:,:,:),4),2))

% all_freq(f_select)*nanmean(pow(:,f_select,1),1)' / sum(nanmean(pow(:,f_select,1),1))

% for ii = 1 : length(i)
%   pp = mp_res(f_select(squeeze(i_res(ii)))-8:f_select(squeeze(i_res(ii)))+8,ii)
%   ff = all_freq(f_select(squeeze(i_res(ii)))-8:f_select(squeeze(i_res(ii)))+8)
%   m_res(ii) = ff*pp / sum(pp);
%   
%   pp = mp_tsk(f_select(squeeze(i_tsk(ii)))-8:f_select(squeeze(i_tsk(ii)))+8,ii)
%   ff = all_freq(f_select(squeeze(i_tsk(ii)))-8:f_select(squeeze(i_tsk(ii)))+8)
%   m_tsk(ii) = ff*pp / sum(pp);
% end
%   
%   m_res = all_freq(f_select(squeeze(i)));
% m_tsk = all_freq(f_select(squeeze(i)));

line([m_res(1) m_res(1)],[min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:))],'color',cmap(1,:),'linewidth',2)
line([m_res(2) m_res(2)],[min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:))],'color',cmap(2,:),'linewidth',2)
line([m_res(3) m_res(3)],[min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:))],'color',cmap(3,:),'linewidth',2)

line([all_freq(f_select(1)) all_freq(f_select(1))],[min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:))],'color','k','linestyle',':')
line([all_freq(f_select(end)) all_freq(f_select(end))],[min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:))],'color','k','linestyle',':')

[~,idx]=min(abs(all_freq(1:4000)-m_tsk(1)));
line([m_tsk(1) m_tsk(1)],[min(mp_tsk(:)) mp_tsk(idx,1)],'color',cmap(4,:),'linewidth',2)

[~,idx]=min(abs(all_freq(1:4000)-m_tsk(2)));
line([m_tsk(2) m_tsk(2)],[min(mp_tsk(:))  mp_tsk(idx,2)],'color',cmap(5,:),'linewidth',2)

[~,idx]=min(abs(all_freq(1:4000)-m_tsk(3)));
line([m_tsk(3) m_tsk(3)],[min(mp_tsk(:)) mp_tsk(idx,3)],'color',cmap(6,:),'linewidth',2)


axis([7 13 min(mp_res(:)) max(mp_res(:))+0.2*max(mp_res(:)) ])

title('Peak frequency (\alpha)');
ylabel('Power'); xlabel('Frequency [Hz]')

subplot(3,2,2)
scatter(fc_all_res(:,1,1),fc_all_res(:,1,2),30,'markerfacecolor',cmap(1,:),'markeredgecolor','w')
lsline; axis square
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([9 11 9 11])

subplot(3,2,4)
scatter(fc_all_res(:,2,1),fc_all_res(:,2,2),30,'markerfacecolor',cmap(2,:),'markeredgecolor','w')
lsline; axis square
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([9 11 9 11])

subplot(3,2,6)
scatter(fc_all_res(:,3,1),fc_all_res(:,3,2),30,'markerfacecolor',cmap(3,:),'markeredgecolor','w')
lsline; axis square
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([9 11 9 11])


subplot(3,2,[5]); hold on
r=corr(fc_all_res(:,:,1),fc_all_res(:,:,2));
imagesc(r,[0.75 1]); axis  square tight
set(gca,'xtick',[1 2 3],'xticklabel',['Pbo1';'Atx1';'Dpz1'])
set(gca,'ytick',[1 2 3],'yticklabel',['Pbo2';'Atx2';'Dpz2'])
h = colorbar; colormap(inferno);
ylabel(h,'Correlation')
% set(gca,'Visible','off')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_power_v%d.pdf',v))

%% PERMUTATION TEST


% pow_res = nanmean(pow_all_res(:,SUBJLIST,:,:),4);
% pow_tsk = nanmean(pow_all_tsk(:,SUBJLIST,:,:),4);

% [~,i] = max(nanmean(pow_res(f_select,:,:),2));
% peak_freq_res = all_freq(f_select(i));
peak_freq_res(1) = tp_peakfreq(nanmean(pow_res(:,:,1),2)',para)
peak_freq_res(2) = tp_peakfreq(nanmean(pow_res(:,:,2),2)',para)
peak_freq_res(3) = tp_peakfreq(nanmean(pow_res(:,:,3),2)',para)

% [~,i] = max(nanmean(pow_tsk(f_select,:,:),2));
% peak_freq_tsk = all_freq(f_select(i));
peak_freq_tsk(1) = tp_peakfreq(nanmean(pow_tsk(:,:,1),2)',para)
peak_freq_tsk(2) = tp_peakfreq(nanmean(pow_tsk(:,:,2),2)',para)
peak_freq_tsk(3) = tp_peakfreq(nanmean(pow_tsk(:,:,3),2)',para)

% -------------------
% PERMUTATION TEST: SETTINGS
% -------------------

para          = [];
para.f        = all_freq;
para.nperm    = 2000;
para.f_select = f_select;
para.all_freq = all_freq;
para.method = 'gaussian'
% -------------------
% PERMUTATION TEST: Atomoxetine vs. placebo (resting state)
% -------------------

dat(:,:,1) = pow_res(:,:,1);
dat(:,:,2) = pow_res(:,:,2);

peak_freq_perm = tp_permtest_peakfreq(dat,para)

perm_diff = peak_freq_perm(:,2)-peak_freq_perm(:,1);
emp_diff = peak_freq_res(:,2)-peak_freq_res(:,1);
  
p(1) = 1-sum(emp_diff<perm_diff)/para.nperm;  
  
% -------------------
% PERMUTATION TEST: Donepezil vs. placebo (resting state)
% -------------------

dat(:,:,1) = pow_res(:,:,1);
dat(:,:,2) = pow_res(:,:,3);

peak_freq_perm = tp_permtest_peakfreq(dat,para)

perm_diff = -diff(peak_freq_perm,[],2);
emp_diff = peak_freq_res(:,3)-peak_freq_res(:,1);
  
p(2) = 1-sum(emp_diff<perm_diff)/para.nperm;  

% -------------------
% PERMUTATION TEST: Task vs rest
% -------------------

dat(:,:,1) = nanmean(pow_res,3);
dat(:,:,2) = nanmean(pow_tsk,3);

peak_freq_perm = tp_permtest_peakfreq(dat,para)

perm_diff = -diff(peak_freq_perm,[],2);
emp_diff = mean(peak_freq_res)-mean(peak_freq_tsk);
  
p(3) = 1-sum(emp_diff<perm_diff)/para.nperm;  

  



%%
isubj = 1;


X = [ones(length(all_freq(200:370)),1) log10(all_freq(200:370)')]
Y = log10(nanmean(dat(200:370,isubj,1),1));

regre = X\Y;

% plot(log10(all_freq(1:380)),log10(all_freq(1:380))*regre(2)+regre(1))
log10(nanmean(dat(:,isubj,1),1))-log10(all_freq(1:380))*regre(2)+regre(1)







