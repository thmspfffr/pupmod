%% pupmod_task_src_power

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
v               = 2;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal_4mm';
f = [8.2500 13.7500; 12.0000 20.0000];
segleng         = 2000;
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v               = 3;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% f = [14 22];
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
end


%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    % %
    if ~exist(sprintf([outdir 'pupmod_task_src_power_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_task_src_power_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    %
    fprintf('Processing s%d m%d f%d ...\n', isubj,m)
    
    for iblock = 1:2
      
      fprintf('Loading MEG data ...\n');
      
      load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      [dat] = megdata2mydata(data); clear data
      
      pars      = [];
      pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      load(pars.sa);
      
      if strcmp(allpara.grid,'cortex_lowres')
        sa.grid_cortex_lowres = select_chans(sa.grid_cortex3000,400);
      end
      
      pars = [];
      
      pars.grid      = allpara.grid;
      
      segleng = 4000;
      fsample = 400;
      all_freq = 0:fsample/segleng:fsample/2;
      
      segshift = segleng / 2;
      maxfreqbin = 400;
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
          pow(iaal,ifoi) = mean(real(diag(cs_src)));
        end
      end
      clear sa
      %
      save(sprintf([outdir 'pupmod_task_src_power_s%d_m%d_b%d_v%d.mat'],isubj,m,iblock,v),'all_freq','pow','f_select');
      
    end
  end
end

error('!')

tp_addpaths
%%
clear fc_all pow_all
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];

ord = pconn_randomization;
peak_range = [7.5 12.5];

v = 2;
load(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],5,1,1,v));


for isubj = SUBJLIST
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 :2
      
      load(sprintf([outdir 'pupmod_src_power_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,v));
      f_select = find(all_freq>=peak_range(1) & all_freq<=peak_range(2));

%       fc_all(isubj,m,iblock) = fc;
      fc_all(isubj,m,iblock) = all_freq(f_select)*mean(pow(:,f_select))' / sum(mean(pow(:,f_select)));
      pow_all(:,isubj,m,iblock) = mean(pow);
      
    end
  end
end

f_select = find(all_freq>=peak_range(1) & all_freq<=peak_range(2));

fc_all  = fc_all(SUBJLIST,:,:)
mean_fc = mean(fc_all,3);
pow_all = mean(pow_all(:,SUBJLIST,:,:),4); 

% fc_all(1) = f*mean(pow_all(:,:,1),2) / sum(mean(pow_all(:,:,1),2));
% fc_all(2) = f*mean(pow_all(:,:,2),2) / sum(mean(pow_all(:,:,2),2));
% fc_all(3) = f*mean(pow_all(:,:,3),2) / sum(mean(pow_all(:,:,3),2));

figure; hold on
subplot(3,2,[1 3]); hold on
co = get(groot,'defaultaxescolororder');
mp = squeeze(mean(pow_all(:,:,:),2));

plot(all_freq(1:400),mp)
axis square 

m = squeeze(mean(mean_fc));

line([m(1) m(1)],[min(mp(:)) max(mp(:))+0.2*max(mp(:))],'color',co(1,:))
line([m(2) m(2)],[min(mp(:)) max(mp(:))+0.2*max(mp(:))],'color',co(2,:))
line([m(3) m(3)],[min(mp(:)) max(mp(:))+0.2*max(mp(:))],'color',co(3,:))

line([all_freq(f_select(1)) all_freq(f_select(1))],[min(mp(:)) max(mp(:))+0.2*max(mp(:))],'color','k','linestyle',':')
line([all_freq(f_select(end)) all_freq(f_select(end))],[min(mp(:)) max(mp(:))+0.2*max(mp(:))],'color','k','linestyle',':')

axis([7 13 min(mp(:)) max(mp(:))+0.2*max(mp(:)) ])

title('Peak frequency (\alpha)');
ylabel('Power'); xlabel('Frequency [Hz]')

subplot(3,2,2)
scatter(fc_all(:,1,1),fc_all(:,1,2),30,'markerfacecolor',co(1,:),'markeredgecolor','w')
lsline; axis square 
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([8.5 10.5 8.5 10.5])

subplot(3,2,4)
scatter(fc_all(:,2,1),fc_all(:,2,2),30,'markerfacecolor',co(2,:),'markeredgecolor','w')
lsline; axis square 
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([8.5 10.5 8.5 10.5])

subplot(3,2,6)
scatter(fc_all(:,3,1),fc_all(:,3,2),30,'markerfacecolor',co(3,:),'markeredgecolor','w')
lsline; axis square 
xlabel('Peak freq #1 [Hz]'); ylabel('Peak freq #2 [Hz]')
axis([8.5 10.5 8.5 10.5])


subplot(3,2,[5]); hold on
r=corr(fc_all(:,:,1),fc_all(:,:,2));
imagesc(r,[0.75 1]); axis  square tight
set(gca,'xtick',[1 2 3],'xticklabel',['Pbo1';'Atx1';'Dpz1'])
set(gca,'ytick',[1 2 3],'yticklabel',['Pbo2';'Atx2';'Dpz2'])
h = colorbar; colormap(inferno); 
ylabel(h,'Correlation')
% set(gca,'Visible','off')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_power_v%d.pdf',v))