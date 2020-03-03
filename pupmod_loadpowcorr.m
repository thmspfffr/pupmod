function fc = pupmod_loadpowcorr(v,SUBJLIST,varargin)

if v == 1
  siz = 91;
elseif v == 12 || v == 19 || v ==23
  siz = 400;
  fc = zeros(siz,siz,length(SUBJLIST),3,2,25,2,'single');
  freqs = [1:25];
elseif v == 20  
  siz = 46;
  fc = zeros(siz,siz,length(SUBJLIST),3,2,25,2,'single');
  freqs = [1:25];
elseif v == 25
  siz = 90;
  freqs = [10 11 12 13];
  fc = zeros(siz,siz,length(SUBJLIST),3,2,length(freqs),2,'single');

else
  error('Invalid version')
end

if varargin{1}==1
  avg = 1;
else
  avg = 0;
end

outdir = '/home/tpfeffer/pupmod/proc/conn/';

addpath ~/pconn/matlab/
ord   = pconn_randomization;
%% LOADING NEW DATA (where on/offset triggers were used)
% -------------------------------------------------------------

ord = pconn_randomization;

clear r_fc
for isubj = 1:length(SUBJLIST)
  fprintf('Loading subject #%d ... \n',isubj)
  
    for m =  1 : 3
      im = find(ord(SUBJLIST(isubj),:)==m);
%       for ifoi = 1 : 13
        
%           load(sprintf('~/pp/proc/conn/pp_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v))
%           fc(:,:,isubj,m,1,ifoi,iblock) = powcorr; clear powcorr
% 
%           load(sprintf('~/pp/proc/conn/pp_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v))
%           fc(:,:,isubj,m,2,ifoi,iblock) = powcorr; clear powcorr
          
%           load(sprintf('~/pp/proc/conn/pp_src_powcorr_test_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
%           fc(:,:,isubj,m,1,:,:) = single(permute(powcorr(:,:,:,freqs),[1 2 4 3])); clear powcorr
% 
%           load(sprintf('~/pp/proc/conn/pp_task_src_powcorr_test_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
%           fc(:,:,isubj,m,2,:,:) = single(permute(powcorr(:,:,:,freqs),[1 2 4 3])); clear powcorr

          load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
          fc(:,:,isubj,m,1,:,:) = permute(powcorr,[1 2 4 3]); clear powcorr

          load(sprintf('~/pupmod/proc/conn/pupmod_task_src_powcorr_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
          fc(:,:,isubj,m,2,:,:) = permute(powcorr,[1 2 4 3]); clear powcorr
        %       end
    
  end
end

if avg == 1
  fc = nanmean(fc(:,:,:,:,:,:,:),7);
else
%   fc = fc(:,:,SUBJLIST,:,:,:,:);
end


