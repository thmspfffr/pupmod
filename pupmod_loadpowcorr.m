function fc = pupmod_loadpowcorr(v,varargin)

if v == 1
  siz = 91;
elseif v == 12
  siz = 400;
elseif v == 20
  siz = 46;
else
  error('Invalid version')
end

if varargin{1}==1
  avg = 1;
else
  avg = 0;
end

outdir = '/home/tpfeffer/pupmod/proc/conn/';

SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/
ord   = pconn_randomization;
%% LOADING NEW DATA (where on/offset triggers were used)
% -------------------------------------------------------------

ord = pconn_randomization;
fc = zeros(siz,siz,34,3,2,13,2,'single');

clear r_fc
for isubj = SUBJLIST
  isubj
  
  for iblock = 1 : 2
    for m =  1 : 3
      im = find(ord(isubj,:)==m);
      for ifoi = 1 : 13
        
          load(sprintf('~/pp/proc/conn/pp_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v))
          fc(:,:,isubj,m,1,ifoi,iblock) = powcorr; clear powcorr

          load(sprintf('~/pp/proc/conn/pp_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v))
          fc(:,:,isubj,m,2,ifoi,iblock) = powcorr; clear powcorr
        
      end
    end
  end
end

if avg == 1
  fc = nanmean(fc(:,:,SUBJLIST,:,:,:,:),7);
else
  fc = fc(:,:,SUBJLIST,:,:,:,:);
end

%% THIS IS THE OLD DATA (where on/offset triggers were not used)
% if avg == 1
% fc = zeros(siz,siz,length(SUBJLIST),3,2,13,'single');
% 
%   for isubj = SUBJLIST
% 
%     fprintf('Processing s%d ... \n',isubj)
% 
%     for ifoi = 1:13
%       for m = 1 : 3
% 
%         im = find(ord(isubj,:)==m);
% 
%         for iblock = 1 : 2
%           clear tmp
%           load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
% 
%           p1(:,:,iblock) = single(powcorr);
% 
%           load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%           p2(:,:,iblock) = single(powcorr);
% 
%         end
% 
%         fc(:,:,isubj,m,1,ifoi) = nanmean(p1,3);
%         fc(:,:,isubj,m,2,ifoi) = nanmean(p2,3);
% 
%         clear p1 p2
%   %       
%       end
%     end
%   end
%   
% else
%   
%   fc = zeros(siz,siz,length(SUBJLIST),3,2,13,2,'single');
% 
%   for isubj = SUBJLIST
% 
%     fprintf('Processing s%d ... \n',isubj)
% 
%     for ifoi = 1:13
%       for m = 1 : 3
% 
%         im = find(ord(isubj,:)==m);
% 
%         for iblock = 1 : 2
%           clear tmp
%           load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
% 
%           fc(:,:,isubj,m,1,ifoi,iblock) = single(powcorr);
%           clear powcorr
%           
%           load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%           fc(:,:,isubj,m,2,ifoi,iblock) = single(powcorr);
%           clear powcorr
%           
%         end
%         
%   %       
%       end
%     end
%   end
%   fc = fc(:,:,SUBJLIST,:,:,:,:);
% 
% end

