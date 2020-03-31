
%% pupmod_src_fcd

% Call pupmod_all_powcorr_periphereal.m next, in order to clean
% estimated FCs from artifacts.

clear
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 23;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'cortex_lowres';
% foi_range       = 2.^[1:.25:7];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% fsample         = 400;
% segleng         = 80;
% segshift        = 40;
% width           = 4;
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
% v               = 24;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% foi_range       = 2.^[1:.25:7];
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 0;
% allpara.tau     = nan;
% weighting       = 1;
% width           = 4;
% --------------------------------------------------------
% VERSION 25 - VOXEL LEVEL, 90 samples AAL (6mm sampling)
% --------------------------------------------------------
v               = 25;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'aal_6mm';
foi_range       = 2.^[1:.25:6];
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
weighting       = 0;
width           = 4;
FOI             = [11];
% --------------------------------------------------------


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

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/
ft_defaults()

outdir   = '/home/tpfeffer/pupmod/proc/fcd/';
addpath /home/tpfeffer/pconn/matlab/

ord   = pconn_randomization;

mask = logical(tril(ones(76,76),-1));

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = 33
  for m = 1
    
%     if ~exist(sprintf([outdir 'pupmod_src_fcd_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' outdir sprintf('pupmod_src_fcd_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end

    for iblock = 2
      
      pars = [];
      pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      sa        = load(pars.sa);
      
      fprintf('Loading MEG data ...\n');
      
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
         continue
      end

     
      if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
        if (data.end_of_recording-600*data.fsample)<1
          data.start_of_recording = 1;
        else
          data.start_of_recording = data.end_of_recording-600*data.fsample;
        end
      elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
        if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
          data.end_of_recording = size(data.trial,2);
        else
          data.end_of_recording = data.start_of_recording+600*data.fsample;
        end
      elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      else
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
      end
      
      data.trial = data.trial(:,data.start_of_recording:data.end_of_recording);
      data.time = data.time(data.start_of_recording:data.end_of_recording);
      
      for ifoi = FOI
        
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)
        data.time(isnan(data.trial(1,:)))=[];
        data.trial(:,isnan(data.trial(1,:)))=[];
        
        cs = tp_wavelet_crossspec(data,foi_range(ifoi),width);
        
        para.iscs = 1;
        para.reg  = 0.05;
        if strcmp(allpara.grid,'aal_4mm')
          sa.sa.L_coarse = sa.sa.L_aal_4mm;
          sa.sa.filt      = pconn_beamformer(real(cs),sa.sa.L_coarse,para);
        else
          sa.sa.filt      = pconn_beamformer(real(cs),sa.sa.L_aal_6mm,para);
        end
        
        [wavelet,~,opt] = tp_mkwavelet(foi_range(ifoi),0.5,400);
        
        % first level segments
        nseg = floor((size(data.trial,2)-opt.n_win)/opt.n_shift+1);
        
        scale=sqrt(nanmean(nanmean(data.trial.^2)));
        data.trial=data.trial/scale;
        
        nanyesno = 0;
        clear datasf1 
        
        for j=1:nseg
          
          dloc2=data.trial(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
          
          if any(isnan(dloc2(:,1)))
            warning('NaN detected')
            datasf1(j,:) = nan(1,size(sa.sa.filt,2));
            nanyesno = 1;
            continue
          end
          
          dataf=dloc2'*wavelet;
          datasf1(j,:)=dataf'*sa.sa.filt;
        end
        
        %2nd level segments
        segleng   = floor(40/(size(wavelet,1)*(1/data.fsample)));
        segshift  = floor(20/(size(wavelet,1)*(1/data.fsample)));
        nseg      = floor((size(datasf1,1)-segleng)/segshift+1);
        
        if nanyesno
          for kk = 1 : size(datasf1,1)
            datasf1(kk,:)=fillmissing(datasf1(kk,:),'spline');
          end
          nanyesno = 0;
        end
        
        for j=1:nseg
          
          fprintf('Segment %d of %d...\n',j,nseg)
          
          dloc2=datasf1((j-1)*segshift+1:(j-1)*segshift+segleng,:);
          
          ns1 = size(datasf1,2);
          res1=zeros(ns1,ns1,'single');
          res2=zeros(ns1,ns1,'single');
          res3=zeros(ns1,ns1,'single');
          res4=zeros(ns1,ns1,'single');
          res5=zeros(ns1,ns1,'single');
          
          for kk = 1 : size(dloc2,1)
            
            dloc = dloc2(kk,:);
            
            for i1=1:ns1
              x1=dloc(i1);
              x2=imag(dloc*conj(x1)./abs(x1));
              y1=abs(x1)^2;
              y2=abs(x2).^2;
              if kk==1
                res1(i1,:)=y1*y2;
                res2(i1,:)=y1;
                res3(i1,:)=y2;
                res4(i1,:)=y1^2;
                res5(i1,:)=y2.^2;
              else
                res1(i1,:)=res1(i1,:)+y1*y2;
                res2(i1,:)=res2(i1,:)+y1;
                res3(i1,:)=res3(i1,:)+y2;
                res4(i1,:)=res4(i1,:)+y1^2;
                res5(i1,:)=res5(i1,:)+y2.^2;
              end
            end
          end
          
          res1=res1/kk;
          res2=res2/kk;
          res3=res3/kk;
          res4=res4/kk;
          res5=res5/kk;
          
          pcorr=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
          pcorr=(pcorr+pcorr')./2;
          if size(pcorr,1)>5000
            reg_idx = 1:90;
            fprintf('Combining atlas regions ...\n')
            for ireg = reg_idx
              ireg
              idx{ireg} = find(sa.sa.aal_label==ireg);
              
            end
            for ireg1 = 1 : size(idx,2)
              for ireg2 = 1 : size(idx,2)
                if ireg1 == ireg2; fc(ireg1,ireg2) = nan; continue; end
                clear tmp
                for iii = 1 : length(idx{ireg1})
                  tmp(iii) = squeeze(mean(pcorr(idx{ireg1}(iii),idx{ireg2}),2));
                end
                fc(ireg1,ireg2) = squeeze(mean(tmp));
              end
            end
            pcorr=fc;
          end
          
          
          para          = [];
          para.transfer = 'to_bcn';
          para.N        = 90;
          
          tmp = tp_match_aal(para,pcorr);
          
          powcorr(:,:,j) = tmp(include_bcn,include_bcn);
          
        end
        
        for iseg = 1 : size(powcorr,3)
          for jseg = 1 : size(powcorr,3)
            tmp1 = powcorr(:,:,iseg);
            tmp2 = powcorr(:,:,jseg);
            fcd(iseg,jseg) = corr(tmp1(mask),tmp2(mask));
          end
        end
        
        clear powcorr
        
        save(sprintf([outdir 'pupmod_src_fcd_s%d_m%d_f%d_b%d_v%d.mat'],isubj,m,ifoi,iblock,v),'fcd');
        
      end
    end
  end
end

error('!')

%%










