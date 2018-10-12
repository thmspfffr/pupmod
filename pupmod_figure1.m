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

