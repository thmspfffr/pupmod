%% % AVALANCHE ANALYSIS
% pupmod_aval

% ----------------------------------------
% VERSION 1
% ----------------------------------------
% v = 1;
% tau = 0.0025:0.0025:0.12;
% T   = 1.5:0.25:5.5;
% ----------------------------------------
% VERSION 2
% ----------------------------------------
% v = 2;
% tau = 0.005:0.0025:0.08;
% T   = 1.5:0.25:5.5;
% ----------------------------------------
% VERSION 3
% ----------------------------------------
v = 3;
tau = 0.005:0.005:0.08;
T   = 1.5:0.25:5.25;
v_preproc = 2;
% ----------------------------------------

% v2: zscore each voxel, v3: zscore across brain

SUBJLIST	= [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
fs        = 400;
dt        = 1/fs;


%%
addpath ~/pconn/matlab/


for isubj = SUBJLIST
  for iblock = 1 : 2
    for m = 1 : 3
      for itau = 1 : length(tau)
        for iT = 1 : length(T)
          
          if ~exist(sprintf(['~/pupmod/proc/pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt'],isubj,iblock,m,itau,iT,v))
            system(['touch ' '~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)]);
          else
            continue
          end
%           
          fprintf('Processing s%d b%d m%d tau%d T%d ...\n',isubj,iblock,m,itau,iT)
          load(sprintf('~/pupmod/proc/pupmod_prepdata_aval_s%d_b%d_m%d_f1_v%d.mat',isubj,iblock,m,v_preproc));
                    
          ns  = round(tau(itau)/dt);
          mm  = mean(mean(dat));
          s   = mean(std(dat));
          z   = zeros(size(dat,1),1);   
          
          nseg = floor(size(dat,1)/ns);

          for ivox = 1 : size(dat,2)
            ivox
            tmp_dat     = [dat(:,ivox)-mm]./s; 
%             dat(:,ivox) = [];
            % zscore each sensor
            tmp         = abs(tmp_dat.*(abs(tmp_dat)>T(iT))); clear tmp_dat
            [~, tmp]    = findpeaks(double(tmp)); 
            z(tmp)      = z(tmp)+1; clear tmp            
          end
          
          clear dat
          
          segshift = ns;
          
          for iseg = 1 : nseg
            
            evts(iseg) = sum(sum(z((iseg-1)*segshift+1:(iseg-1)*segshift+ns,:)));
            
          end
          
          clear z
           
          a = 0; br = 0; clear aval
          
          fprintf('Identifying avalanches ... \n')
          
          while any(evts>0)
            
            start = find(evts~=0,1,'first');
           
            a = a + 1;
            
            aval(a) = evts(start);
            
            evts(start)=0;
            
            if (start+1)>=length(evts)-1
              break
            end
            
            i = 1;
            
            while evts(start+i)>0
              
              aval(a) = aval(a) + evts(start+i);
              
              evts(start+i)=0;
              
              i = i + 1;      
              
              if (start+i)==length(evts)
                break
              end
              
            end
            
            par.dur(a) = i;
            
            i = 0;
            
          end
          
          fprintf('Identifying avalanches ... Done. \n')

          [par.a,par.b]=hist(aval,1:20:100000);
          
          save(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,m,itau,iT,v)],'par');
          
          clear aval evts par evts aval a z
          
        end
      end
    end
  end
end


error('!')

%% CLEAN FILES
v = 3;
tau = 0.005:0.005:0.08;
T   = 1.5:0.25:5.25;
v_preproc = 2;

outdir   = '/home/tpfeffer/pconn/proc/dfa/';

cnt = 0;
v = 3;
cnt_exist = 0;
SUBJLIST  = [31 32 33 34];
for m = 1 : 3
  for isubj = SUBJLIST
    for itau = 1:length(tau)
      for iT = 1 : length(T)
        for iblock = 1 : 2
        
      
          if exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,m,itau,iT,v)]) && exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)])
            cnt_exist = cnt_exist + 1;

            continue
          elseif exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,m,itau,iT,v)]) && ~exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)])
            system(['touch ' '~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)]);

          elseif exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)]) && ~exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,m,itau,iT,v)])
            warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
            delete(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)])
            cnt = cnt + 1;
          elseif ~exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)]) && exist(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,m,itau,iT,v)])
            system(['touch ' '~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d_processing.txt',isubj,iblock,m,itau,iT,v)]);
          else
            warning('Nothing exists')
            cnt = cnt+1;
          end
      
        end
      end
    end
  end
end
cnt

%%
clear slp

ord       = pconn_randomization;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v         = 3;

for isubj =  SUBJLIST
  isubj
  for iblock = 1 : 2
    for m = 1 : 3
      im = find(ord(isubj,:)==m);
      
      
      for itau = 1 : length(tau)
        for iT = 1 : length(T)
          
          load(['~/pupmod/proc/' sprintf('pupmod_aval_s%d_b%d_m%d_tau%d_T%d_v%d.mat',isubj,iblock,im,itau,iT,v)]);

          tmp = par.a;
          idx = find(tmp~=0);
          
          s=pconn_regress(log10(par.b(idx(2:end))),log10(tmp(idx(2:end)))');
          
          if any(isnan(s))
            error('!')
          end
          slp(isubj,iblock,m,itau,iT) = s(2);
        end
      end
    end
  end
end

slp = squeeze(nanmean(slp(SUBJLIST,:,:,:,:),2));

%%

m_slp = squeeze(nanmean(slp(:,:,:,:,:),1));
% all_slp = squeeze(nanmean(slp(:,:,:,:,:),2));

for i = 1 : 16
  for j = 1 : 16
    
    [~,p(i,j)]=ttest(slp(:,2,i,j)-slp(:,1,i,j));
  end
end


%% COMPUTE CORELATIONS WITH DFA
v = 2;
cnt = 0;
for isubj  = SUBJLIST
  cnt = cnt + 1;
  for m = 1  : 3
    for ifoi = 1 : 4
      
      im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v));
      
      d(isubj,m,ifoi) = nanmean(par.dfa(:));
      
      
    end
  end
end

d = d(SUBJLIST,:,:);
%%
ss = squeeze(nanmean(slp,2));
dd = squeeze(nanmean(d,2));

ifoi = 3;
for i = 1 : size(slp,3)
  for j = 1 : size(slp,4)
    
    [r(i,j) p(i,j)] = corr(squeeze(-ss(:,i,j)),squeeze(dd(:,ifoi)));
    
  end
end

r(isnan(r)) = 0;
r(isnan(p)) = 0;

%
figure_white;

k(1)=subplot(2,3,1)
imagesc(flipud((r.*(p<0.05))'),[-1 1]); colormap(jet)
xlabel('tau [ms]'); ylabel('T [SD]');
title('Corr(beta,alpha)');
axis square

par1 = -squeeze(nanmean(slp(:,2,:,:)))';
par2 = -squeeze(nanmean(slp(:,1,:,:)))';
t    = squeeze(ttest(slp(:,2,:,:),slp(:,1,:,:),'dim',1,'alpha',0.08))
t(isnan(t)) = 0;

k(2)=subplot(2,3,2)
imagesc(flipud(squeeze(t)').*(flipud(par1)-flipud(par2)),[-0.05 0.05]); colormap(jet)

xlabel('tau [ms]'); ylabel('T [SD]');
title('Atx vs. pbo');
axis square

k(3)=subplot(2,3,5)
imagesc((flipud(par1)-flipud(par2)),[-0.05 0.05]); colormap(jet)
xlabel('tau [ms]'); ylabel('T [SD]');
title('Atx vs. pbo');
axis square

par1 = -squeeze(nanmean(slp(:,3,:,:)))';
par2 = -squeeze(nanmean(slp(:,1,:,:)))';
[t]    = squeeze(ttest(squeeze(slp(:,3,:,:)),squeeze(slp(:,1,:,:)),'dim',1,'alpha',0.08));
t(isnan(t)) = 0;

k(4)=subplot(2,3,3)
imagesc(flipud(squeeze(t)').*(flipud(par1)-flipud(par2)),[-0.05 0.05]); colormap(jet)

xlabel('tau [ms]'); ylabel('T [SD]');
title('Dpz vs. pbo');
axis square

k(5)=subplot(2,3,6)
imagesc((flipud(par1)-flipud(par2)),[-0.05 0.05]); colormap(jet)
xlabel('tau [ms]'); ylabel('T [SD]');
title('Dpz vs. pbo');
axis square
% slp=squeeze(slp);
% 
% figure_white;
% 
% imagesc(flipud(-slp'),[0 2]); colormap(jet)
% xlabel('tau [ms]'); ylabel('T [SD]');
% 
para.cond = 'bth'

behav = nanmean(pconn_getbehavior(para),2);

for i = 1 : size(slp,3)
  for j = 1 : size(slp,4)
    
    [r(i,j) p(i,j)] = corr(squeeze(-ss(:,i,j)),behav);
    
  end
end

r(isnan(r)) = 0;

k(6)=subplot(2,3,4)
imagesc(flipud(r'),[-0.5 0.5]); colormap(jet)
xlabel('tau [ms]'); ylabel('T [SD]');
title('Corr(behav,alpha)');
axis square

set(k,'tickdir','out','xtick',3:8:length(tau),'xticklabel',tau(3:8:end).*1000);
set(k,'tickdir','out','ytick',1:4:length(T),'yticklabel',T(end:-4:1));

print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_aval_f%d_v%d.jpeg',ifoi,v));
%%

s = squeeze(nanmean(slp,2));
% % 

