%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pupmod_prepdata_aval

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v               = 1;
% v_grid          = 2; % 4 = aal
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% allpara.filt    = 'eloreta';
% allpara.grid    = 'xcoarse';
% v_rawdata       = 6;
% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v               = 2;
v_grid          = 3; % 4 = aal
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
allpara.filt    = 'eloreta';
allpara.grid    = 'cortex';
v_rawdata       = 6;
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
outdir   = '/home/tpfeffer/pupmod/proc/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for ifoi = 1 : 1
  for m = 1 : 3
    for isubj = SUBJLIST
%       
      if ~exist(sprintf([outdir 'pupmod_prepdata_aval_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pupmod_prepdata_aval_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      
      for iblock = 1:2
        
        disp(sprintf('Loading MEG data ...'));

        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
        
        [dat,epleng] = megdata2mydata(data_low); clear data_low
        
        pars = [];
        
        pars.fsample   = 400;
        pars.segleng   = 400;
        pars.segshift  = 400;
        pars.epleng    = size(dat,1);
        
        cs = data2cs_event(dat,pars.segleng,pars.segshift,pars.epleng,40,pars.fsample);
        
        % get spatial filter
        pars = [];
        pars.grid = allpara.grid;
        pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        pars.filt = allpara.filt;
        pars.cs   = cs;
        pars.foi  = [1 40];
        
        % get spatial filter as defined elsewhere
        filt      = get_spatfilt(pars);
               
        dat = dat*filt;
                    
        save(sprintf([outdir 'pupmod_prepdata_aval_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'dat','-v7.3');
        
        clear par
        
      end
      
    end
  end
end

error('STOP')


%% % AVALANCHE ANALYSIS
fs  = 400;
dt  = 1/fs;
tau = 0.02;


ns = tau/dt;

load(sprintf(['~/pupmod/proc/dat4etienne_meg_s%d_b%d_m%d_f%d_v%d.mat'],4,1,1,1,2));

z=zeros(size(dat));

for ivox = 1 : size(dat,2)
  
  tmp       = abs(dat(:,ivox).*(abs(zscore(dat(:,ivox)))>3));
  [~, tmp]  = findpeaks(tmp);
  z(tmp,ivox) = 1;
  
end

nseg = floor(size(dat,1)/ns);
segshift = ns;

for iseg = 1 : nseg
  
  evts(iseg) = sum(sum(z((iseg-1)*segshift+1:(iseg-1)*segshift+ns,:)));

end

a = 0; br = 0;

while any(evts>0)
    
  start = find(evts~=0,1,'first');
  
  a = a + 1;
  
  aval(a) = evts(start);
  
  evts(start)=0;
  
  i = 1;
  
  while evts(start+i)>0 & ~br
    
    aval(a) = aval(a) + evts(start+i);
    
    evts(start)=0;
    
    i = i + 1;
    
    if (start+i)>length(evts)
      br = 1;
      break
    end
    
  end
  
  i = 0;
    
end
    
  
  
  

  
  




