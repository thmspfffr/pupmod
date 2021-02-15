function fc = pupmod_loadpowcorr(v,SUBJLIST,varargin)


if varargin{1}==1
  avg = 1;
else
  avg = 0;
end

if v == 1 || v==2 || v == 3 || v == 4 | v == 333
  siz = 400;
  if avg == 0
    fc = zeros(siz,siz,length(SUBJLIST),3,2,17,2,'single');
  else
    fc = zeros(siz,siz,length(SUBJLIST),3,2,17,'single');
  end
elseif v == 33
  siz = 90;
  freqs = [4 5 6 7 8 9 10];
  if avg == 0
    fc = zeros(siz,siz,length(SUBJLIST),3,2,length(freqs),2,'single');
  else
    fc = zeros(siz,siz,length(SUBJLIST),3,2,length(freqs),1,'single');
  end
else
  error('Invalid version')
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
  
  for m =  1 :3
    im = find(ord(SUBJLIST(isubj),:)==m);
    
    load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
    
    if avg == 1
      fc(:,:,isubj,m,1,:) = tanh(nanmean(atanh(permute(powcorr.*1.73,[1 2 4 3])),4)); clear powcorr
    else
      fc(:,:,isubj,m,1,:,:) = permute(powcorr.*1.73,[1 2 4 3]); clear powcorr
    end
    %
    load(sprintf('~/pupmod/proc/conn/pupmod_task_src_powcorr_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
    %
    if avg == 1
      fc(:,:,isubj,m,2,:) = tanh(nanmean(atanh(permute(powcorr.*1.73,[1 2 4 3])),4)); clear powcorr
    else
      fc(:,:,isubj,m,2,:,:) = permute(powcorr.*1.73,[1 2 4 3]); clear powcorr
    end
    
  end
end
