load /home/tpfeffer/pconn/proc/preproc/pconn_postproc_s9_m3_b1_v6.mat
data = data_low.trial{1}(1,:);
%%

figure; set(gcf,'color','w'); hold on

plot(zscore(data(10000:10400)),'k','linewidth',2); axis tight off

print(gcf,'-dpdf',sprintf('~/pupmod_figure1a.pdf'))


 flp = 9;           % lowpass frequency of filter
 fhi = 13;
 
 delt  = 1/400;            % sampling interval
 k     = 4;                  % 2nd order butterworth filter
 fnq   = 1/(2*delt);       % Nyquist frequency
 Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
 [bfilt,afilt] = butter(k,Wn);
 
 % COMPUTE DFA, EXTRACT HURST
 % ---------------------------
 env =filtfilt(bfilt,afilt,data);
 figure; set(gcf,'color','w'); hold on
plot(zscore(data(10000:10400)),'color',[0.8 0.8 0.8],'linewidth',1); axis tight off

 plot(zscore(env(1,10000:10400)),'r','linewidth',2); axis tight off
 
print(gcf,'-dpdf',sprintf('~/pupmod_figure1b.pdf'))

env2 = abs(hilbert(env));
figure; set(gcf,'color','w'); hold on

 plot(env2(1,10000:10400),'m','linewidth',2); axis tight off
  plot(env(1,10000:10400),'r','linewidth',2); axis tight off

print(gcf,'-dpdf',sprintf('~/pupmod_figure1c.pdf'))

%% CORRELATION WITH SC
% load connectome
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));
C = C(include_bcn,include_bcn);

mask = logical(tril(ones(size(C,1),size(C,1)),-1));
foi_range       = unique(round(2.^[1:.5:7]))

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',1));

for ifoi = 1 : 13
  
  tmp  = nanmean(cleandat(include_bcn,include_bcn,:,1,1,ifoi),3);
  fc = tmp(mask);
  [r(ifoi) p(ifoi)] = corr(fc,C(mask));
  
end

figure; set(gcf,'color','w')
subplot(2,1,1)
area(r,'linestyle','none')
set(gca,'XTick',1:13,'xticklabel',foi_range)
xlabel('Frequency [Hz]')
ylabel('Correlation with SC')
axis([1 13 0 0.4])
% axis square

print(gcf,'-dpdf',sprintf('~/pupmod_correlationwithSC.pdf'))

