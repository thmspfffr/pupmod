%% COMPUTE POWER IN AAL ATLAS
% pupmod_src_pow
% performs distance-weighting as implemented for FC data
% see Brookes et al. (2016) NeuroImage

% last updated 25-09-2018

clear
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% -----------------------------------------------------
% VERSiON 2
% -----------------------------------------------------
v       = 1;
FOI     = [11 16];
smoo    = [2.75 4];
% -----------------------------------------------------

tp_addpaths

indir   = '/home/tpfeffer/pupmod/proc/src/';
outdir   = '/home/tpfeffer/pupmod/proc/src/';

%%

for ifoi = 1: length(FOI)
  for m = 1 : 3
    for isubj = SUBJLIST

      fn = sprintf('pupmod_src_pow_s%d_m%d_f%d_v%d',isubj,m,ifoi,v);
      if tp_parallel(fn,'~/pupmod/proc/src/',1)
        continue
      end

      clear d
      
      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
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
            out.pow(:,2) = nan(90,1); 
          else
            out.pow(:,1) = nan(90,1); 
          end
        end
        
        fprintf('Processing MEG s%dm%df%d%b...\n',isubj,m,ifoi,iblock);

        if length(d) == 2
          load(['~/pconn/proc/preproc/' d(iblock).name]);
        else
          load(['~/pconn/proc/preproc/' d(1).name]);
        end
        
        cfg = [];
        cfg.length =  1/0.2;
        cfg.overlap = 0;
        
        data      = ft_redefinetrial(cfg,data);
        
        cfg               = [];
        cfg.method        = 'mtmfft';
        cfg.output        = 'powandcsd';
        cfg.taper         = 'dpss';
        cfg.pad           = 'nextpow2';
        cfg.foi           = mean(FOI(ifoi,:));
        cfg.keeptrials    = 'no';
        cfg.tapsmofrq     = smoo(ifoi);
        
        [~,  csd]         = ft_freqanalysis(cfg, data);
   
        clear data
        
        % watch out! L_coarse is correct even for cortex grid
        load(['~/pconn/proc/src/' sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,7)]);        
        pos = sa.grid_aal4mm_indi;
        
        for ireg = 1 : 90
          
          idx = sa.aal_label==ireg;

          lf = squeeze(sa.L_aal_4mm(:,idx,:));
          para = [];
          para.reg = 0.05;
          para.iscs = 1;
          filt 	= pconn_beamformer(csd,lf,para);
          
          if length(idx)>1
            com = mean(pos(idx,:));
          else
            com = pos(idx,:);
          end

          dist = sqrt((com(1)-pos(idx,1)).^2 + (com(2)-pos(idx,2)).^2 + (com(3)-pos(idx,3)).^2);

          % In cm or mm?
          xx = pos; for i=1:3; xx(:,1)=xx(:,i)-mean(xx(:,i)); end
          xrad = mean(sqrt(sum(xx.^2,2)));

          if xrad > 20
            w    = exp((-dist.^2)./400);
          else
            w    = exp((-10*dist.^2)./400);
          end
          
          % computed weighted power values
          out.pow(ireg,iblock)         = sum(w.*diag(filt'*real(csd)*filt))/length(idx);
          
        end

        clear csd filt
        
      end
      
      save(sprintf([outdir 'pupmod_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'out');
      
      tp_parallel(fn,'~/pmod/proc/ei/',0)
      
      clear out
      
    end
  end
end

error('STOP')

%% CLEAN NON PROCESSED FILES
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap  = cbrewer('div', 'RdBu', 500,'pchip');
% cmap  = cmap(end:-1:1,:);
% TEST PLOT
load sa_meg_template;

grid  = sa_meg_template.grid_cortex3000;
g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = sa_meg_template.grid_cortex3000;
% par_interp = spatfiltergauss(ans,g1,dd,g2);
%%
para = [] ;
para.colorlimits = [min(par_interp1) max(par_interp1)];

% PLOT RESULTS
tp_showsource(par_interp1,cmap,sa_meg_template,para);

%% PLOT ALL
ifoi = 4
v = 2

ord	= pconn_randomization;
cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    
      im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn/proc/src/' 'pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
      
      pow_all1(:,isubj,m) = nanmean(par.pow,2); clear par
      
    	load(sprintf(['~/pconn/proc/src/' 'pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,2));
      pow_all2(:,isubj,m) = nanmean(par.pow,2);

      
  end
end

pow1 = pow_all1(:,SUBJLIST,:);
pow2 = pow_all2(:,SUBJLIST,:);


% d=;
% for isubj = 1 : size(pow,2);
%   d = nanmean(nanmean(pow(:,isubj,:),2),3);


%%
icontr = 2;
contrast = [2 1; 3 1; 2 3]; 

% d = (nanmean(pow(:,:,3),2)-nanmean(pow(:,:,1),2))./(nanmean(pow(:,:,2),2)+nanmean(pow(:,:,1),2));
[t2,p2,~,s1]=ttest(pow1(:,:,contrast(icontr,1)),pow1(:,:,contrast(icontr,2)),'dim',2);
d1 = s1.tstat; %d = d.*t2;

[t1,p1,~,s2]=ttest(pow2(:,:,contrast(icontr,1)),pow2(:,:,contrast(icontr,2)),'dim',2);
d2 = s2.tstat; %d = d.*t2;

para.nperm = 10000;
para.tail = 0;
para.paired=1;
para.neigh = get_neighbours(g1);
para.method = 'dependentT'
para.alpha = 0.025;
para.clusteralpha = 0.05;
para.minneigh = 2;

% s = tp_clusterperm(pow(:,:,[2 1]),para)

  par_interp = spatfiltergauss(d1,g1,dd,g2);
  
  para = [] ;
  para.colorlimits = [-1.96 1.96];
%   para.colorlimits = [min(par_interp) max(par_interp)];

  % PLOT RESULTS
  tp_showsource(par_interp,cmap,sa_meg_template,para); drawnow
  
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pow_c%d_f%d_v%d.mat.jpg',icontr,ifoi,1))

par_interp = spatfiltergauss(d2,g1,dd,g2);
  
  para = [] ;
  para.colorlimits = [-1.96 1.96];
%   para.colorlimits = [min(par_interp) max(par_interp)];

  % PLOT RESULTS
  tp_showsource(par_interp,cmap,sa_meg_template,para); drawnow
  
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pow_c%d_f%d_v%d.mat.jpg',icontr,ifoi,2))


  %%
for isubj = 1 : 28
    p1 = pow1(:,isubj,2)-pow1(:,isubj,1);
    p2 = pow2(:,isubj,2)-pow2(:,isubj,1);
  
    r(isubj)=corr(p1,p2)
%     par_interp = spatfiltergauss(p,g1,dd,g2);
%     para.colorlimits = [log10(min(par_interp)) log10(max(par_interp))];
%     tp_showsource(log10(par_interp),hot,sa_meg_template,para); drawnow
% print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pow_s%d_f%d_v%d.mat.jpg',isubj,ifoi,v))
end
  
  
  
% end


