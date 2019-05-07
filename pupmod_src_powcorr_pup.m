%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pupmod_src_powcorr_pup

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v               = 1;
v_grid          = 4; % 4 = aal
fsample         = 400;
SUBJLIST        = [4 5 6 7 9 10 11 12 13 15 16 19 20 21 22 23 24];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal';
foi_range       = unique(round(2.^[1:.5:7]));
para.smo        = foi_range./4;
para.segleng    = 1 ./ para.smo;
para.epleng     = 60;
lpc             = 0;
v_pup           = 10;
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v               = 2;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.smo        = foi_range./4;
% para.segleng    = 1 ./ para.smo;
% para.epleng     = 60;
% lpc             = 0;
% v_pup           = 10;
% para.wavelet    = 'ft';
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v               = 3;
% v_grid          = 4; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.smo        = foi_range./4;
% para.segleng    = 1 ./ para.smo;
% para.epleng     = 60;
% lpc             = 0;
% v_pup           = 10;
% para.wavelet    = 'bp_filt';
% para.scnd_filt = 0;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
addpath ~/pconn/matlab

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir  = '/home/tpfeffer/pupmod/proc/';
freq    = 1;

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for ifoi = 1 : length(foi_range)
  for m = 1 : 3
    for isubj = SUBJLIST
      
      if ~exist(sprintf([outdir 'pupmod_src_powcorr_pup_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_src_powcorr_pup_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      %
      % ---------------------------------------------------
      % READ IN PUPIL
      % ---------------------------------------------------
      
      d=dir(['/home/tpfeffer/pconn/proc/pup/' sprintf('pconn_postproc_pupdfa_s%d_b*_m%d_v%d.mat',isubj,m,v_pup)]);
      
      if v_pup>9
        a=11;
      else
        a=10;
      end
      
      if length(d)<1
        num_blocks(isubj,m) = 0;
%         continue
      elseif length(d)==1 && str2double(d(1).name(end-a)) == 1
        num_blocks(isubj,m) = 1;
        bl = 1;
        warning(sprintf('only one block which is b%d',bl));
      elseif  length(d)==1 && str2double(d(1).name(end-a)) == 2
        num_blocks(isubj,m) = 2;
        bl = 2;
      else
        num_blocks(isubj,m) = 2;
        bl = 1;
      end
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      
      for iblock = bl : num_blocks(isubj,m)
        
        disp(sprintf('Loading MEG data ...'));

        load(sprintf('~/pconn/proc/pup/pconn_postproc_pupdfa_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_pup))

        if size(dat,2)>size(dat,1)
          dat = dat';
        end
        
        pup.dil=tp_interp_blinks(pup.dil);
                
        pars = [];
        
        pars.fsample   = 400;
        pars.segleng   = round(para.segleng(ifoi).*fsample);
        pars.segshift  = round(fsample*para.segleng(ifoi)/2);
        pars.foi       = foi_range(ifoi);
        pars.epleng    = size(dat,1);
        
        cs = data2cs_wavelet(dat,pars.segleng,pars.segshift,pars.epleng,pars.foi,pars.fsample);
        
        % get spatial filter
        pars = [];
        pars.grid = allpara.grid;
        pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        pars.filt = allpara.filt;
        pars.cs   = cs;
        pars.foi  = [foi_range(ifoi) foi_range(ifoi)];
        
        get spatial filter as defined elsewhere
        filt      = get_spatfilt(pars);               
        sa        = load(pars.sa);
        
        pars = [];
        
        pars.fsample   = 400;
        pars.segleng   = round(para.segleng(ifoi).*fsample);
        pars.segshift  = round(fsample*para.segleng(ifoi)/2);
        pars.foi       = foi_range(ifoi);
        
        pars.epleng    = para.epleng*pars.fsample;
        pars.epshift   = pars.epleng-round(0.95*pars.epleng);
       	pars.grid      = 'medium';
        pars.aal       = 0;
        pars.wavelet   = para.wavelet;
        pars.scnd_filt = para.scnd_filt; 
        % COMPUTE CONNECTIVITY OVER TIME
        nep = floor((size(dat,1)-pars.epleng)/pars.epshift+1);

        for iep = 1 : nep

          timerange             = [(iep-1)*pars.epshift+1 (iep-1)*pars.epshift+pars.epleng];
          pup_dat               = pup.dil(timerange(1):timerange(2));
          tmp_dat               = dat(timerange(1):timerange(2),:);   
          
          
          if ~lpc
            tmp_powcorr= tp_powcorr_ortho(tmp_dat,pars,sa);
          else
            tmp_powcorr = tp_data2lpc_jackknife(tmp_dat,pars,filt,filt);
          end
          
          par.powcorr(:,:,iep) = tp_match_aal(pars,tmp_powcorr);
          
          clear tmp_powcorr
        
        % compute various things here:
        % kura order, pairwise kura, power corr, pupil, 

          par.pupmean(iep)      = mean(pup_dat);
          par.pupstd(iep)     	= std(pup_dat);
          par.pup(:,iep)       	= pup_dat;

          
        end
        
        % COMPUTE MEAN CONNECTIVITY
        pars.epleng    = size(dat,1);
        powcorr = tp_powcorr_ortho(dat,pars,sa);
        powcorr = tp_match_aal(pars,powcorr);
        pow = nonzeros(triu(powcorr,1));
      
        for iep = 1 : nep
          for jep = 1 : nep 
            
            par.fcd(iep,jep) = corr(nonzeros(triu(par.powcorr(:,:,iep),1)),nonzeros(triu(par.powcorr(:,:,jep),1)));
          
          end
        end

        par.fcdvar = var(nonzeros(triu(par.fcd,1)));
                
        save(sprintf([outdir 'pupmod_src_powcorr_pup_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'par','-v7.3');
        
        clear par mydata dat
        
      end    
    end
  end
end

error('STOP')


%%
clear fcd
v = 1;
figure_white;
for isubj = SUBJLIST
  for iblock = 1 : 1
    for m = 1 : 1
      foicnt = 0;
      for ifoi = 1:13
        
        foicnt = foicnt + 1;
        load(sprintf([outdir 'pupmod_src_powcorr_pup_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v));
%         
        
        

% subplot(3,4,foicnt)
%         imagesc(par.fcd,[0 .5]); colormap(jet)
%         axis square
%         
%         fcd(:,:,foicnt) = par.fcd;
        
      end
      
    end
  end
end
%%


ifoi = 7;

for it = 1 : length(fcd)
  [r(it),p(it)]=corr(fcd(:,it,ifoi),par.pupmean');
end

t = 1 : 50
starting_lambda = 2;


figure_white; hold on;

plot(find(p<0.05),zeros(sum(p<0.05),1),'k.')
plot(r,'linewidth',3)

%% AUTOCORR DECAY

ifoi = 7
for i=1:100%length(fcd)
  
  k(:,i)=acf(fcd(i:end,i,ifoi),50);
  f = @(lambda) sum((k(:,i)'-exp(-lambda.*t)).^2); % das ist die funktion. innerhalb der function ist die exp-function und drumherum der summed squared error (sse). diese funktion wird miminiert
  [l_hat(i),fval] = fminsearch(f,[starting_lambda]);
  
end

figure_white; hold on;
plot(zscore(par.pupmean(1:100)));
plot(zscore(l_hat))

%% FC DIFFERENCES
fois = [1 3 4 5 7 9 13];
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% figure_white;
ord = pconn_randomization;

% powcorr = zeros(

for ifoi = 1 : 13
  clear fcd powcorr
  v = 1;
  cnt = 0;
  for isubj = SUBJLIST
    isubj
      for m = 1 : 3
        for iblock = 1 : 2
 
        im = find(ord(isubj,:)==m);

        foicnt = 0;        

%         cnt = cnt + 1
        try
          load(sprintf([outdir 'pupmod_src_powcorr_pup_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,im,ifoi,v));
        catch me
          continue
        end
        
        idx = par.pupmean>median(par.pupmean);
        powcorr(:,:,1,isubj,m,iblock)=mean(par.powcorr(:,:,idx),3);
        idx = par.pupmean<median(par.pupmean);
        powcorr(:,:,2,isubj,m,iblock)=mean(par.powcorr(:,:,idx),3);
        
        
      end
    end
  end


  p = nanmean(powcorr(:,:,:,SUBJLIST,:,:),6);
  
  p1 = squeeze(nanmean(p(:,:,1,:,:),4));
  p2 = squeeze(nanmean(p(:,:,2,:,:),4));
  
  s11 = std(nonzeros(triu(p1(:,:,1))));
  s12 = std(nonzeros(triu(p1(:,:,2))));

  s21 = std(nonzeros(triu(p1(:,:,1))));
  s22 = std(nonzeros(triu(p1(:,:,3))));
  
  figure_white; title(sprintf('Freq: %d',ifoi));
  imagesc(nanmean(p1,3)-nanmean(p2,3),[-0.01 0.01]); colormap(jet)
  title(sprintf('Freq: %d',ifoi));
 
  figure_white; 
  
  subplot(2,4,1)
  imagesc(p1(:,:,2)-p1(:,:,1),[-0.02 0.02]); colormap(jet);axis square
  subplot(2,4,2)
  imagesc(p2(:,:,2)-p2(:,:,1),[-0.02 0.02]); colormap(jet);axis square
  subplot(2,4,3)
  imagesc([p1(:,:,2)-p1(:,:,1)] - [p2(:,:,2)-p2(:,:,1)],[-0.02 0.02]); colormap(jet);axis square
  
  subplot(2,4,5)
  imagesc(p1(:,:,3)-p1(:,:,1),[-0.02 0.02]); colormap(jet);axis square
  subplot(2,4,6)
  imagesc(p2(:,:,3)-p2(:,:,1),[-0.02 0.02]); colormap(jet);axis square
  subplot(2,4,7)
  imagesc([p1(:,:,3)-p1(:,:,1)] - [p2(:,:,3)-p2(:,:,1)],[-0.02 0.02]); colormap(jet);axis square
  title(sprintf('Freq: %d',ifoi));

  
end
   



        


% figure_white; hold on;



