%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_pow

% requires execution of pconn_src_compute_common_filter.m first.
% this is where the common dipole orientation is determined.

clear
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% -----------------------------------------------------
% VERSiON 12 - cortex
% -----------------------------------------------------
v       = 12;
foi_range       = unique(round(2.^[1:.5:7]));
bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
% -----------------------------------------------------

tp_addpaths

indir   = '/home/tpfeffer/pupmod/proc/pow/';
outdir   = '/home/tpfeffer/pupmod/proc/pow/';

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid, 'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'genemaps')
  v_grid = 13;
elseif strcmp(allpara.grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(allpara.grid,'cortex800')
  v_grid = 16;
end
%%

for ifoi = 6
  for m = 1 : 3
    for isubj = SUBJLIST
      
%       if ~exist(sprintf([outdir 'pupmod_src_pow_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
%         system(['touch ' outdir sprintf('pupmod_src_pow_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
%       else
%         continue
%       end
      
      clear d
      
      fprintf('Processing s%d m%d f%d ..\n', isubj,m,ifoi)
      
      d=dir(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
      
      if length(d)==1
        if v < 10
          blocks = str2num(d(1).name(end-7));
        else
          blocks = str2num(d(1).name(end-8));
        end
      elseif length(d) == 2
        blocks = [1 2];
      end
            
      for iblock = blocks
        
        if length(blocks) == 1
          if iblock == 1
            par.pow(:,2) = nan(400,1); 
          else
            par.pow(:,1) = nan(400,1); 
          end
        end
        
        fprintf('Processing MEG s%dm%df%d%b...\\n',isubj,m,ifoi,iblock);
            
        if length(d) == 2
          load(['~/pconn/proc/preproc/' d(iblock).name]);
        else
          load(['~/pconn/proc/preproc/' d(1).name]);
        end
        
%         [dat] = megdata2mydata(data); clear data
%         
%         cfg = [];
%         cfg.length =  1/0.2;
%         cfg.overlap = 0;
%         
%         data      = ft_redefinetrial(cfg,data);
        
%         % COMPUTE CROSS SPECTRAL DENSITY, using multitaper method
%         % ------------------------------
%         cfg               = [];
%         cfg.method        = 'mtmconvol';
%         cfg.output        = 'powandcsd';
%         cfg.taper         = 'hanning';
%         cfg.pad           = 'nextpow2';
%         cfg.foi           = foi_range(ifoi);
%         cfg.keeptrials    = 'no';
%         cfg.tapsmofrq     = (bpfreq(ifoi,2)-bpfreq(ifoi,1))/2;
%         
%         [~,  csd]         = ft_freqanalysis(cfg, data);  
        % ------------------------------
        data = megdata2mydata(data); 
        segleng = 1600;
        segshift = 1600;
        epleng = size(data,1);
        maxfreqbin = 300;
        [cs] = data2cs_event(data,segleng,segshift,epleng,maxfreqbin,para);
        
        f = 0:400/1600:200;
        
        csd = mean(real(cs(:,:,f>=bpfreq(6,1) & f<=bpfreq(6,2))),3);
%         clear data
        
        % LOAD HEAD/FORWARD MODEL
        % ------------------------------
        pars      = [];
        pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        sa        = load(pars.sa);
        % ------------------------------
        [A, Ao, po]=mkfilt_lcmv(sa.sa.L_coarse,csd,.05);
        
        noisecov = reshape(sa.sa.L_coarse,[size(sa.sa.L_coarse,1) size(sa.sa.L_coarse,2)*3])*reshape(sa.sa.L_coarse,[size(sa.sa.L_coarse,1) size(sa.sa.L_coarse,2)*3])';
        
        para        = [];
        para.iscs   = 1;
        para.reg    = 0.05;
%         para.noisecov = noisecov;
        [filt,pow]  = pconn_beamformer(csd,sa.sa.L_coarse,para);
    
        save(sprintf([outdir 'pupmod_src_pow_s%d_m%d_f%d_b%d_v%d.mat'],isubj,m,ifoi,iblock,v),'pow');

        clear csd filt pow 
        
      end
      
      
      clear par
      
    end
  end
end

error('STOP')

%% PLOT
ord   = pconn_randomization;
for ifoi = 6
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for isubj = SUBJLIST
      for iblock = 1 : 2
        load(sprintf([outdir 'pupmod_src_pow_s%d_m%d_f%d_b%d_v%d.mat'],isubj,im,ifoi,iblock,v));
        
        allpow(:,isubj,m,iblock) = pow;
      end
    end
  end
end

allpow= allpow(:,SUBJLIST,:,:)
      
  %%    
      
load /home/gnolte/meth/templates/sa_template.mat;

pow = log10(nai)
par = pow;
para = [];
para.clim = [min(pow) max(pow)];
para.cmap = hot;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/test.png');
tp_plot_surface(pow,para)


