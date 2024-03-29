%% pupmod_src_conn

clear all

% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v         = 2;
% fsample   = 400;
% v_postproc = 6;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% v_grid    = 1; % 1 = coarse; 2 = coarse; 3 = cortex
% allpara.filt = 'eloreta';
% allpara.grid = 'coarse';
% foi_range = unique(round(2.^[1:.5:7]));
% i_calc = [1 100];
% i_fit 	= [2 100];
% para.smo = foi_range./4;
% para.segleng = 1 ./ para.smo;
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
v         = 3;
v_cs      = 2; % cross spec: 2 = wavelets
v_postproc = 6;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
v_grid    = 3; % 1 = coarse; 2 = coarse; 3 = cortex
allpara.filt = 'jh_lcmv';
allpara.grid = 'cortex';
foi_range = unique(round(2.^[1:.5:7]));
i_calc = [1 100];
i_fit 	= [2 100];
para.smo = foi_range./4;
para.segleng = 1 ./ para.smo;
% --------------------------------------------------------

% foi_range = foi_range(foi_range<40);

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pupmod/proc/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 1;

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    for ifoi = 1:length(foi_range)
      
      if ~exist(sprintf([outdir 'pupmod_src_conn_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_conn_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      
      for iblock = 1:2
        
        disp(sprintf('Loading MEG data ...'));
        
        load(['~/pconn/proc/src/' sprintf('pconn_sens_cs_s%d_m%d_b%d_f%d_v%d.mat',isubj,m,iblock,ifoi,v_cs)]);
        
        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_postproc));
        data_low.trial{1} = data_low.trial{1} + data_hi.trial{1}; clear data_hi
        
        [dat,epleng] = megdata2mydata(data_low); clear data_low
        dat = dat';

        % get spatial filter
        pars = [];
        pars.grid = allpara.grid;
        pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        pars.filt = allpara.filt;
        pars.cs   = cs;
        pars.foi  = [foi_range(ifoi) foi_range(ifoi)];
        % get spatial filter as defined elsewhere
        filt      = get_spatfilt(pars);
        
        pars = [];
        
        pars.fsample   = 400;
        pars.segleng   = round(para.segleng(ifoi).*fsample);
        pars.segshift  = round(fsample*para.segleng(ifoi)/2);
        pars.foi       = foi_range(ifoi);
        pars.epleng    = 8000;
        pars.aal       = 0;
        
        % COMPUTE POWER CORRELATIONS
        % dat should be n x nchans
        tmp = tp_orthopowercorr(dat',pars,filt,filt);
          
          % Convert to AAL
        powcorr = single(tp_grid2aal(tmp,pars));              
        
        save(sprintf([outdir 'pupmod_src_conn_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        
      end
    end
  end
end
  
  
  error('!')
  
  
%%

FOI = [2 4; 4 8; 8 13; 13 38; 54 128];

addpath ~/pconn/matlab/
  
ord = pconn_randomization;

for ifoi = 1 : 9
  
%   foi = [find(foi_range==FOI(ifoi,1)):find(foi_range==FOI(ifoi,2))];

  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      for iblock = 1 : 2
        
        load(sprintf([outdir 'pupmod_src_conn_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));

        s(:,:,isubj,m,ifoi,iblock) = mean(par.powcorr,3);

      end
    end
  end
end


%%
figure; set(gcf,'color','white');



clim = [0 0.25];

s_plt = nanmean(s.*1.77,6);
s_plt = s_plt(:,:,SUBJLIST,:,:);

ss = squeeze(nanmean(nanmean(nanmean(nanmean(s_plt,4),3),2),1));

%%











































%% PLOT PEAK FREQ OF VAR
a = squeeze(nanmean(s_plt,3));
a = squeeze(nanmean(a,3));
a = squeeze(nanmean(a,2));

for i = 1 : 91
%   for j = 1 : 91
    [~,t(i)]=max(a(i,:));
    idx(i) = foi_range(t(i));
%   end
end

%% PLOT ON AAL COARSE

load sa_meg_template;
mri   = sa_meg_template.mri;

addpath ~/pconn/matlab/
load aalmask_grid_cortex3000
grid = sa_meg_template.grid_cortex3000;

dat = zeros(3000,1);

for i = 1 : 91
  
  aalidx = find(aalgrid.mask==i);
  dat(aalidx) = idx(i);
  
end


g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

vc = sa_meg_template.vc;

viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

dd = .1;
%%
par_interp = spatfiltergauss(dat,g1,dd,g2);
   
figure; set(gcf,'color','white'); hold on;

para = [] ;
% r = max(abs(min(d)),abs(max(d)));
%   para.colorlimits = [min(d(d>eps)) max(d)];
para.colorlimits = [0 38];

for iplot = 1 : 6
  
  subplot(3,2,iplot)
  
  para.myviewdir = viewdir(iplot,:);
  a = sa_meg_template.cortex10K;
  
  if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
  elseif iplot == 6
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
  end
  
  pconn_showsurface(a,para,par_interp)
  
  colormap(parula)
  
  camlight headlight
  
end

% saveas(gcf,sprintf([plotdir 'pconn_src_contrast_f%d_v%d.fig'],ifoi,v),'fig')

%%













































%%

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)

ifoi = 6;


clim = [0 0.25];

s_plt = squeeze(nanmean(nanmean(nanmean(s(:,:,SUBJLIST,:,ifoi,:),6),5),3));

dat = zeros(3000,3);

sss = squeeze(nanmean(s_plt,2));

for i = 1 : 91
  
  aalidx = find(aalgrid.mask==i);
  dat(aalidx,:) = repmat(sss(i,:),[length(aalidx) 1]);
  
end

dat = dat(:,2)-dat(:,1);

g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

vc = sa_meg_template.vc;

viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

dd = .1;

par_interp = spatfiltergauss(dat,g1,dd,g2);
   
figure; set(gcf,'color','white'); hold on;

para = [] ;
% r = max(abs(min(d)),abs(max(d)));
%   para.colorlimits = [min(d(d>eps)) max(d)];
para.colorlimits = [-max([abs(min(dat)) abs(max(dat))]) max([abs(min(dat)) abs(max(dat))])];

for iplot = 1 : 6
  
  subplot(3,2,iplot)
  
  para.myviewdir = viewdir(iplot,:);
  a = sa_meg_template.cortex10K;
  
  if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
  elseif iplot == 6
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
  end
  
  pconn_showsurface(a,para,par_interp)
  
  colormap(cmap)
  
  camlight headlight
  
end



%%

for iblock = 1 : 2
  
  load(sprintf('/home/tpfeffer/pupmod/proc/pupmod_src_powcorr_s4_m1_b%d_f10_v1.mat',iblock))


  for i = 1  : 91
    for j = 1 : 91

      [r(i,j,iblock),p(i,j,iblock)]=corr(squeeze(par.powcorr(i,j,:)),par.pup');

    end
  end

end




  
  
  
  %% GET TIME COURSES IN ORDER TO COMPUTE LOCAL BIFURCATION PARAMETER
  % SOURCE SETTINGS
  para.grid = 'coarse';
  para.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
  para.filt = 'lcmv';
  para.fs   = 400;
  para.foi  = [8 12];
  para.cs   = cs;
  para.aal  = 0;
  
  [dat A1]  = pupmod_src_timecourse(mydata,para);
  
  
  
  
  
  
  
  cd ~/Downloads/
  
  load SC_AAL_90.mat
  
  areas = 90;
  C = SC/max(max(SC))*0.2;
  
  time_steps = 600;
  
  a = 0.2;
  we = 0.8;
  sig = 0.01;
  omega = 10;
  
  threshold = 0;
  polarity = 1;
  
  dt = 0.0025;
  
  xs = zeros(time_steps/2,areas);
  x = 0.1 * randn(areas,1);
  y = 0.1 * randn(areas,1);
  nn = 0;
  
  Isubdiag = find(tril(ones(areas),-1));
  
  
  for t=0:dt:500
    t
    sumax = C*x-diag(C*repmat(x',areas,1));
    sumay = C*y-diag(C*repmat(y',areas,1));
    
    x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
    y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
    
  end
  
  for t=0:dt:5000  %32000
    t
    % if t > 500
    %     a = 0.1;
    % end
    
    sumax = C*x-diag(C*repmat(x',areas,1));
    sumay = C*y-diag(C*repmat(y',areas,1));
    
    x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
    y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
    
    
    if mod(t,2)==0
      nn=nn+1;
      xs(nn,:)=x';
      ri=sqrt(x.*x+y.*y);
      rimax=max(ri);
      ri=ri/rimax;
      kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
    end
    
  end
  
  
  FC_simul = corrcoef(xs(1:nn,:));
  
  
  % RUN DFA
  f_win = [1 100];
  c_win = [0.5 150];
  
  siginfo = nbt_Info;
  siginfo.converted_sample_frequency =  5;
  
  t = nbt_doDFA(abs(hilbert(xs)), siginfo, [1 100],[0.5 150],0.5,0,0,[]);
  % cc = corrcoef(FC_emp(Isubdiag),FC_simul(Isubdiag));
  v = var(abs(hilbert(xs)),[],1);
  m = mean(abs(hilbert(xs)));
  
  
  %%
  % -0.2
  m = 0.0117
  v = 3.78e-05
  d = 0.56
  
  % 0
  m = 0.271
  v = 2.08e-04
  d = 0.92
  
  % 0.2
  m = 0.44;
  v = 8.72e-05;
  d = 0.56
  
  save('~/pmod/proc/pmod_hopf_par2.mat','t','v','m')
  
  
  Corr = cc(1,2);
  
  Corr_sim_mean = mean(mean(FC_simul));
  
  data = xs';
  [Peak,Amp] = Threshold_fMRI(data, areas, threshold, polarity);
  
  figure;
  subplot(2,1,1)
  pcolor(Peak)
  subplot(2,1,2)
  plot(1:time_steps/2+1, xs(:,1:66))
  xlim([0 time_steps/2+1])
  
