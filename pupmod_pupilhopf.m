%% pupmod_pupilhopf

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v          = 6;
fsample    = 400;
para.bif   = 0;
para.G     = 1.5;
para.sig   = .01;
para.freq  = 10;
para.ts    = 100; % timesteps
plt = 0;
nruns = 10000;
% --------------------------------------------------------

outdir   = '/home/tpfeffer/pupmod/proc/';


%% RUN HOPF MODEL
if ~exist(sprintf('~/pupmod/proc/pupmod_pupilhopf_noise_v%d.mat',v))
  noise = rand(para.ts/(1/fsample)+1,nruns/2);
  noise = 2*(noise-repmat(mean(noise),[para.ts/(1/fsample)+1 1]));
  save(sprintf('~/pupmod/proc/pupmod_pupilhopf_noise_v%d.mat',v),'noise');
else
  load(sprintf('~/pupmod/proc/pupmod_pupilhopf_noise_v%d.mat',v));
end
  
for irun = 1 : nruns
  
  if ~exist(sprintf('~/pupmod/proc/pupmod_pupilhopf_r%d_v%d_proc.txt',irun,v))
    system(['touch ' sprintf('~/pupmod/proc/pupmod_pupilhopf_r%d_v%d_proc.txt',irun,v)]);
  else 
    continue
  end

  fprintf('Processing run %d ...\n', irun);
  
  if irun < nruns/2+1
    % dynamic modulation
    gmod = noise(:,irun);
  else
%   no modulation
    gmod = ones(para.ts/(1/fsample)+1,1).*mean(noise(:,irun-nruns/2));
  end
  
  for isig = 1 : length(para.sig)
    for ibi = 1 : length(para.bif)
      for iG = 1 : length(para.G)
        
        fs = 400;
        
        load ~/Downloads/SC_AAL_90.mat
        
        areas       = 90;
        C           = SC/max(max(SC))*0.2;
        
        time_steps  = para.ts;
        
        a           = para.bif(ibi);  % bifurcation param.
        we          = para.G(iG);     % coupling
        sig         = para.sig(isig);           % amp. of noise?
        f           = para.freq;             % freq. = omega
        
        dt          = 1/fsample;
        
        amp         = zeros(time_steps/dt,areas);
        ys          = zeros(time_steps/dt,areas);
        xs          = zeros(time_steps/dt,areas);
        
        kura        = zeros(time_steps/dt,1);
        
        x           = 0.1 * randn(areas,1);
        y           = 0.1 * randn(areas,1);
        nn          = 0;
        omega       = 2*pi*f;
        
        Isubdiag = find(tril(ones(areas),-1));
        
        fprintf('Initializing...\n');
        
        xs = zeros(para.ts/(1/fsample)+1,90);
        ys = zeros(para.ts/(1/fsample)+1,90);
        amp = zeros(para.ts/(1/fsample)+1,90);
        kura = zeros(para.ts/(1/fsample)+1,1);
        
        for t=0:dt:50
          
          sumax = C*x-diag(C*repmat(x',areas,1));
          sumay = C*y-diag(C*repmat(y',areas,1));
          
          x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
          y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
          
        end
        
        if plt
          figure; set(gcf,'color','w');
        end
        fprintf('Running model...\n');
        ts_cnt = 0;
        for t=0:dt:time_steps  %32000
          ts_cnt = ts_cnt+1;
%           fprintf('Time step %.3f / %d ...\n',t,time_steps);
          
          sumax = C*x-diag(C*repmat(x',areas,1));
          sumay = C*y-diag(C*repmat(y',areas,1));
          
          %       x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
          %       y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
          x = x+dt*((gmod(ts_cnt)+a).*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
          y = y+dt*((gmod(ts_cnt)+a).*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
          
          %       if mod(t,2)==0
          nn=nn+1;
          xs(nn,:)=x';
          ys(nn,:)=y';
          
          amp(nn,:) = sqrt(x.*x+y.*y);
          %       if plt
          %         if nn >  100
          %           subplot(2,1,1)
          %           plot(xs(nn:-1:nn-100),ys(nn:-1:nn-100),'-')
          %           axis square tight
          %           %       axis([-sqrt(para.bif)-sqrt(para.bif)*0.2 sqrt(para.bif)+sqrt(para.bif)*0.2 -sqrt(para.bif)-sqrt(para.bif)*0.2 sqrt(para.bif)+sqrt(para.bif)*0.2]); grid on;
          %           subplot(2,1,2)
          %           plot(101:-1:1,amp(nn:-1:nn-100,1))
          %           axis square tight
          %           drawnow
          %         end
          %       end
          
          ri=amp(nn,:)';
          rimax=max(ri);
          ri=ri/rimax;
          kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
          
        end
        
        fc_sim(:,:) = corrcoef(amp(1:nn,:));
        
        nchan = 90;
        f = 11;
        seglen  = 400; % in samples
        nseg    = floor(size(xs,1)/seglen);
        overlap = 0.5;
        win = hanning(seglen);
        cs=zeros(nchan,nchan,nseg);
        
        for iseg = 1 : nseg
          
          idx = (iseg-1)*seglen*overlap+1:(iseg-1)*seglen*overlap+seglen;
          
          dat = xs(idx,:);
          
          datfft = fft(repmat(win,[1 size(dat,2)]).*dat);
          
          cs(:,:,iseg)      = cs(:,:,iseg)+conj(datfft(f,:)'*datfft(f,:));
          %   abscoh(:,:,iseg)  = cs(:,:,iseg)./sqrt(diag(cs(:,:,iseg))*diag(cs(:,:,iseg))');
          
          
        end
        
        cs = nanmean(cs,3);
        
        coh = [abs(cs).^2]./[abs(diag(cs))*abs(diag(cs))'];
        save(sprintf('~/pupmod/proc/pupmod_pupilhopf_r%d_v%d.mat',irun,v),'coh','fc_sim')

      end
    end
  end
end


error('!')


%%
clear m_coh
clear m_amp
nruns = 10000; v = 6;

for irun = 1 : nruns
irun
  load(sprintf('~/pupmod/proc/pupmod_pupilhopf_r%d_v%d.mat',irun,v));
  
  m_coh(irun) = mean(coh(:));
  m_amp(irun) = mean(fc_sim(:));
end

m1 = mean(m_coh(1:nruns/2));
m2 = mean(m_coh(nruns/2+1:end));

s1 = std(m_coh(1:nruns/2))/sqrt(nruns/2);
s2 = std(m_coh(nruns/2+1:end))/sqrt(nruns/2);

figure; set(gcf,'color','w');

plot(1:2,[m1 m2],'k.','markersize',20);

line([1 1],[m1-s1 m1+s1],'color','k','linewidth',0.5);
line([2 2],[m2-s2 m2+s2],'color','k','linewidth',0.5);

% axis([0 3 .2 .3])

[~,p]=ttest(m_coh(1:nruns/2),m_coh(nruns/2+1:end))

m1 = mean(m_amp(1:nruns/2));
m2 = mean(m_amp(nruns/2+1:end));

s1 = std(m_amp(1:nruns/2))/sqrt(nruns/2);
s2 = std(m_amp(nruns/2+1:end))/sqrt(nruns/2);

figure; set(gcf,'color','w');

plot(1:2,[m1 m2],'k.','markersize',20);

line([1 1],[m1-s1 m1+s1],'color','k','linewidth',0.5);
line([2 2],[m2-s2 m2+s2],'color','k','linewidth',0.5);

% axis([0 3 .090 .105])

[~,p]=ttest(m_amp(1:nruns/2),m_amp(nruns/2+1:end))

%%

struct(:,:,:,1) = coh;
struct(:,:,:,2) = fc_sim;
struct(:,:,:,3) = repmat(SC,[1 1 21]);

save('~/ampcorr_coh_sc.mat','struct');

%%


figure; set(gcf,'color','white'); hold on

for i = 1 : 45
  
  a = fc_sim(:,:,i);
  b = coh(:,:,i);
  
  r(i) = corr(a(:),b(:));
  
  ra(i) = corr(a(:),SC(:));
  rb(i) = corr(b(:),SC(:));
  
  m(i) = nanmean(a(:));
  n(i) = nanmean(b(:));
  
end

figure; set(gcf,'color','w'); hold on


plot(para.sig(1:45),ra,'r','linewidth',3)
plot(para.sig(1:45),rb,'b','linewidth',3)

xlabel('Noise');
ylabel('Correlation(amp. corr,abs(coh))');


figure; set(gcf,'color','w'); hold on

plot(para.sig(1:45),r,'b','linewidth',3)

xlabel('Noise');
ylabel('Correlation(amp. corr,abs(coh))');

figure; set(gcf,'color','w'); hold on
plot(m,'r','linewidth',3)
plot(n,'b','linewidth',3)

xlabel('Noise');
ylabel('Mean connectivity')


%%
figure; set(gcf,'color','white');

% PLOT CONNECTIVITY MATRICES

subplot(2,3,1);

imagesc(coh);  axis square;

subplot(2,3,2);

imagesc(fc_sim);  axis square;

subplot(2,3,3);

icoh = imag([(cs).^2]./[diag(cs)*diag(cs)']);

imagesc(icoh);  axis square;

% PLOT CORRELATIONS

subplot(2,3,4);

scatter(coh(:),fc_sim(:)); axis square;  box on;
r=corr(coh(:),fc_sim(:));

subplot(2,3,5);

scatter(coh(:),icoh(:)); axis square;  box on;
r=corr(coh(:),icoh(:));

subplot(2,3,6);

scatter(fc_sim(:),icoh(:)); axis square; box on;
r=corr(fc_sim(:),icoh(:));


fprintf('Amp. corr with SC: %.2f\n',corr(fc_sim(:),SC(:)));

fprintf('Coh. with SC: %.2f\n',corr(coh(:),SC(:)));


