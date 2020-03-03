function fc = pupmod_loadfcd(v,SUBJLIST,varargin)

if v == 1
%   siz = 91;
% elseif v == 12 || v == 19 || v ==23
%   siz = 400;
%   fc = zeros(siz,siz,length(SUBJLIST),3,2,25,2,'single');
%   freqs = [1:25];
% elseif v == 20  
%   siz = 46;
%   fc = zeros(siz,siz,length(SUBJLIST),3,2,25,2,'single');
%   freqs = [1:25];
elseif v == 25
  freq = [ 11  ];
%   siz = 90;
%   fc = zeros(siz,siz,length(SUBJLIST),3,2,3,2,'single');
%   freqs = [10 11 12];
else
%   error('Invalid version')
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
fc = nan(100,100,length(SUBJLIST),3,1);
ord = pconn_randomization;

clear r_fc
for isubj = 1:length(SUBJLIST)
  isubj
  
  for iblock = 1 : 2
    for m =  1 : 3
      im = find(ord(SUBJLIST(isubj),:)==m);
      for ifoi = 1:1
        ifreq = freq(ifoi);
        for iblock = 1 : 2
          
          try
          load(sprintf('~/pupmod/proc/fcd/pupmod_src_fcd_s%d_m%d_f%d_b%d_v%d.mat',SUBJLIST(isubj),im,ifreq,iblock,v))
          fc(1:size(fcd,1),1:size(fcd,1),isubj,m,1) = single(fcd); clear fcd
          catch me
            continue
          end
%           size(fcd)
          try
          load(sprintf('~/pupmod/proc/fcd/pupmod_src_task_fcd_s%d_m%d_f%d_b%d_v%d.mat',SUBJLIST(isubj),im,ifreq,iblock,v))
          fc(1:size(fcd,1),1:size(fcd,1),isubj,m,2) = single(fcd); clear fcd
          catch me
            continue
          end
%           size(fcd)
        end
      end
    end
  end
end

% if avg == 1
%   fc = nanmean(fc(:,:,:,:,:,:,:),7);
% else
% %   fc = fc(:,:,SUBJLIST,:,:,:,:);
% end

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


