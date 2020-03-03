% pupmod_sens_peakfreq


clear
% ---------
% VERSION 1: gaussian fit
% ---------
v = 1
method = 'gaussian';
detrend = 0;
% ---------
% VERSION 2: center of mass
% ---------
% v = 2;
% method = 'com';
% detrend = 0;
% ---------
% VERSION 3: center of mass
% ---------
% v = 3;
% method = 'gaussian';
% detrend = 1;
% ---------

addpath ~/Documents/MATLAB/fieldtrip-20160919/; ft_defaults

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
outdir = '~/pupmod/proc/sens/';
%%
for isubj = SUBJLIST
  for m = 1:3
    
    if ~exist(sprintf([outdir 'pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_sens_peakfreq_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    for iblock =1:2
      
      fprintf('Processing s%d m%d b%d ...\n', isubj,m,iblock)

      try
        load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        continue
      end
      out.cnt =0;
      % TEST FILTER
      % --------------------------
%       all_idx = find(isnan(data.trial(1,:)));
%       start = 1;
%       while 1
% %         
%         if start == 1
%           idx=find(isnan(data.trial(1,start:end)),1,'first');
%         else
%           idx=all_idx(find(all_idx>start,1,'first'));
%           if isempty(idx)
%             idx=size(data.trial,2);
%           end
%         end
%         if isnan(data.trial(1,idx))
%           tmp_dat = data.trial(:,start:idx-1);
%         else
%           tmp_dat = data.trial(:,start:idx);
%         end
%         
%         if size(tmp_dat,2)<300
%           warning('data short')
%           data.trial(:,start:idx-1) = nan;
%           out.cnt = out.cnt+1;
%           while isnan(data.trial(1,idx))
%             idx=idx+1;
%             if idx==size(data.trial,2)
%               break
%             end
%           end
%           
%           if idx==size(data.trial,2)
%             break
%           end
%           start=idx;
%           continue
%         end
%         
%         clear mirr_dat
%         
%         mirr_dat=padarray(tmp_dat',800,'symmetric','both');       
%         mirr_dat_filt = ft_preproc_highpassfilter(mirr_dat',400,2,4,'but');
%         tmp_dat = mirr_dat_filt(:,801:size(mirr_dat_filt,2)-800);
%         if isnan(data.trial(1,idx))
%           data.trial(:,start:idx-1)=tmp_dat;
%         else
%           data.trial(:,start:idx)=tmp_dat;
%         end
%         
%         while isnan(data.trial(1,idx))
%           idx=idx+1;
%           if idx==size(data.trial,2)
%             break
%           end
%         end
%         
%         if idx==size(data.trial,2)
%           break
%         end
%         start=idx;
%       end
%      % --------------------------
      
      dat = data.trial';
      
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
      
      dat = dat(data.start_of_recording:data.end_of_recording,:);
      dat(isnan(dat(:,1)),:)=[];
      
      
%       for ichan = 1 : size(dat,2)
%                   ichan
%         dx = abs(fft(dat(:,ichan))).^2; dx = dx(1:length(dx)/2+1);
%         para.f = 0:400/size(dat,1):400/2;
%         para.f_select = para.f >7 & para.f <13;
%         para.method = method;
%         para.detrend = 0;
%         para.detrendfreq = [1 100];
%         dx_resample = resample(dx,10,400);
%         out.f = para.f(1:40:end);
%         out.pow7_13hz(:,ichan) = dx_resample;
%         [out.pf(ichan)] = tp_peakfreq(dx,para);
%         
%       end
      
%       dx = mean(abs(fft(dat)).^2,2);
%       [out.pf_avg out.fun_avg] = tp_peakfreq(dx,para);
      
      [px,f]=pwelch(dat,hanning(4000),0.5,0:0.1:50,400,'power');
      out.pwelch = px;
      out.pwelchfreq = f;
      
      save(sprintf('~/pupmod/proc/sens/pupmod_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v),'out');
      
      clear out dx_resample out.f
      
    end
  end
end

% COMPUTE REST

for isubj = SUBJLIST
  for m = 1:3
    
    if ~exist(sprintf([outdir 'pupmod_rest_sens_peakfreq_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pupmod_rest_sens_peakfreq_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    
    for iblock = 1:2
      
       fprintf('Processing s%d m%d b%d ...\n', isubj,m,iblock)

      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        out.cnt = 0;
        %           out.pf =
        %           save(sprintf([outdir 'pupmod_sens_peakfreq_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        continue
      end
      
      % TEST FILTER
      % --------------------------
%       all_idx = find(isnan(data.trial(1,:)));
%       start = 1;
%       while 1
% %         
%         if start == 1
%           idx=find(isnan(data.trial(1,start:end)),1,'first');
%         else
%           idx=all_idx(find(all_idx>start,1,'first'));
%           if isempty(idx)
%             idx=size(data.trial,2);
%           end
%         end
%         if isnan(data.trial(1,idx))
%           tmp_dat = data.trial(:,start:idx-1);
%         else
%           tmp_dat = data.trial(:,start:idx);
%         end
%         
%         if size(tmp_dat,2)<300
%           warning('data short')
%           data.trial(:,start:idx-1) = nan;
%           out.cnt = out.cnt+1;
%           while isnan(data.trial(1,idx))
%             idx=idx+1;
%             if idx==size(data.trial,2)
%               break
%             end
%           end
%           start=idx;
%           if idx==size(data.trial,2)
%             break
%           end
%           continue
%         end
%         
%         mirr_dat=padarray(tmp_dat',800,'symmetric','both');       
%         mirr_dat_filt = ft_preproc_highpassfilter(mirr_dat',400,2,4,'but');
%         tmp_dat = mirr_dat_filt(:,801:size(mirr_dat_filt,2)-800);
% 
%          if isnan(data.trial(1,idx))
%           data.trial(:,start:idx-1)=tmp_dat;
%         else
%           data.trial(:,start:idx)=tmp_dat;
%         end
% 
%         while isnan(data.trial(1,idx))
%           idx=idx+1;
%           if idx==size(data.trial,2)
%             break
%           end
%         end
%        
%         if idx==size(data.trial,2)
%           break
%         end
%          start=idx;
%       end
      % --------------------------
      
      dat = data.trial';
      
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
      
      dat = dat(data.start_of_recording:data.end_of_recording,:);
      dat(isnan(dat(:,1)),:)=[];
      
%       for ichan = 1 : size(dat,2)
%                   ichan
%         dx = abs(fft(dat(:,ichan))).^2; dx = dx(1:length(dx)/2+1);
%         para.f = 0:400/size(dat,1):400/2;
%         para.f_select = para.f >7 & para.f <13;
%         para.method = 'gaussian';
%         para.detrend = 0;
%         para.detrendfreq = [1 100];
%         dx_resample = resample(dx,10,400);
%         out.f = para.f(1:40:end);
%         out.pow7_13hz(:,ichan) = dx_resample;
%         [out.pf(ichan)] = tp_peakfreq(dx,para);
%         
%       end
      
%       dx = mean(abs(fft(dat)).^2,2);
%       [out.pf_avg out.fun_avg] = tp_peakfreq(dx,para);
      
       [px,f]=pwelch(dat,hanning(4000),0.5,0:0.1:50,400,'power');
      out.pwelch = px;
      out.pwelchfreq = f;
      save(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v),'out');
      
      clear out
      
    end
  end
end


%%
addpath ~/pconn/matlab
clear pf_all pf_avg_all pow
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
clear ps 
v = 2
ord = pconn_randomization;
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      im = find(ord(isubj,:)==m);
      try
        load(sprintf('~/pupmod/proc/sens/pupmod_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
%         
%         ii = 0; cnt=1; clear ff
%         while 1        
%           if isempty(find(out.f>0+ii&out.f<0.2+ii))
%             break
%           elseif (out.f>0+ii)>=200
%             break
%           end
%           ff(cnt)=mean(mean(out.pow7_13hz(out.f>0+ii&out.f<0.2+ii,:),2));
%           cnt = cnt+1;
%           ii=ii+0.2;
%         end

          pow(:,isubj,m,2,iblock)=mean(out.pwelch,2);
%         if size(out.pf,2)==274
%           pf_all(:,isubj,m,2,iblock) = out.pf;
%         else
%           pf_all(:,isubj,m,2,iblock) = pconn_sens_interp274_alena(isubj,out.pf,im);
%         end
%         pf_avg_all(isubj,m,2,iblock) = out.pf_avg;
%          [~,i]=max(smoothdata(mean(out.pow7_13hz(para.f>8&para.f<12),1),'gaussian', 30));
%         k=para.f(para.f>8&para.f<12);
%         pf_new(isubj,m,2,iblock) = k(i);
%         pow(:,isubj,m,2,iblock)=ff;
      catch me
%         pf_all(:,isubj,m,2,iblock) = nan(274,1);
%         pf_avg_all(isubj,m,2,iblock)  = nan;
%         pf_new(isubj,m,2,iblock) =nan;
%         pow(isubj,m,2,iblock)=nan;
pow(:,isubj,m,2,iblock)=nan(501,1);
      end

%       end
    end
    
    for iblock = 1 : 2
      try
        load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_peakfreq_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v));
        
%         ii = 0; cnt=1; clear ff
%         while 1        
%           if isempty(find(out.f>0+ii&out.f<0.2+ii))
%             break
%           elseif (out.f>0+ii)>=200
%             break
%           end
%           ff(cnt)=mean(mean(out.pow7_13hz(out.f>0+ii&out.f<0.2+ii,:),2));
%           cnt = cnt+1;
%           ii=ii+0.2;
%         end
pow(:,isubj,m,1,iblock)=mean(out.pwelch,2);
%         pow(:,isubj,m,1,iblock)=ff;
%         if size(out.pf,2)==274
%           pf_all(:,isubj,m,1,iblock) = out.pf;
%           ps(:,isubj,m,1,iblock) = mean(out.pow7_13hz,2);
%         else
%           pf_all(:,isubj,m,1,iblock) = pconn_sens_interp274_alena(isubj,out.pf,im);
%           ps(:,isubj,m,1,iblock) = pconn_sens_interp274_alena(isubj,mean(out.pow7_13hz,2),im);
%         ps(:,isubj,m,1,iblock) = mean(out.pow7_13hz,2);
%         end
%         pf_avg_all(isubj,m,1,iblock) = out.pf_avg;
%          [~,i]=max(smoothdata(mean(out.pow7_13hz(para.f>8&para.f<12),1),'gaussian', 30));
%         k=para.f(para.f>8&para.f<12);
%         pf_new(isubj,m,1,iblock) = k(i);


%         for 
%         pow(:,isubj,m,1,iblock)=nanmean(out.pow7_13hz(1:1500,:),2);
      catch me
%         pf_all(:,isubj,m,2,iblock) = nan(274,1);
%         pf_avg_all(isubj,m,2,iblock)  = nan;
%         pf_new(isubj,m,2,iblock) =nan;
%         pow(isubj,m,2,iblock)=nan;
pow(:,isubj,m,1,iblock)=nan(501,1);
      end
      clear out
    end
  end
end

pf_all = nanmean(pf_all(:,SUBJLIST,:,:,:),5);
pf_avg_all= nanmean(pf_avg_all(SUBJLIST,:,:,:),4);

%% PLOT AVERAGE POWER SPECTRUM (ACROSS SUBJECTS)
smo = 20;
figure; set(gcf,'color','w'); hold on

plot(f,smooth(nanmean(pow(:,:,1,1),2),smo))
plot(f,smooth(nanmean(pow(:,:,1,2),2),smo))

idx=out.pwelchfreq>7&out.pwelchfreq<14;
[~,i1]=max(smooth(mean(pow(idx,:,1,1),2),smo));
[~,i2]=max(smooth(mean(pow(idx,:,1,2),2),smo));
ff=out.pwelchfreq(idx);
pf_diff=ff(i1)-ff(i2);

nperm=10000;
all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);

for iperm = 1 : nperm
  
    idx1 = all_idx1(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permpow(:,i,1) = pow(:,i,1,idx1(i));
      permpow(:,i,2) = pow(:,i,1,idx2(i));     
    end
    
    idx=out.pwelchfreq>7&out.pwelchfreq<14;
    [~,i1]=max(smooth(mean(permpow(idx,:,1),2),smo));
    [~,i2]=max(smooth(mean(permpow(idx,:,2),2),smo));
    ff=out.pwelchfreq(idx);
    
    perm_pf(iperm,1)=ff(i1);
    perm_pf(iperm,2)=ff(i2);
    
end


perm_diff= perm_pf(:,1)-perm_pf(:,2);

p=sum((pf(1)-pf(2))< perm_diff)/nperm

%%
para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_cnt = behav;

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);

%% correlate drug effects
ipharm = 2;
iblock = 1:2;
ifoi = 6;

tmp_fc = nanmean(nanmean(fc(:,:,:,:,2,ifoi,iblock),7),6);

d_beh = nanmean(behav_cnt(ipharm,:,iblock)-behav_cnt(1,:,iblock),3);
d_pf = nanmean(b(:,ipharm,iblock)-b(:,1,iblock),3);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7)-nanmean(nanmean(nanmean(nanmean(fc(:,:,:,1,2,ifoi,iblock),1),2),7),6));
% 
% d_beh = (nanmean(behav_cnt(ipharm,:,iblock)-behav_cnt(1,:,iblock),3))./nanmean(behav_cnt(1,:,iblock),3);
% d_pf = (nanmean(b(:,ipharm,iblock)-b(:,1,iblock),3))./nanmean(b(:,1,iblock),3);
% d_meg = (nanmean(reshape(tmp_fc(:,:,:,ipharm),[400*400 28]))-nanmean(reshape(tmp_fc(:,:,:,1),[400*400 28])))./nanmean(reshape(tmp_fc(:,:,:,1),[400*400 28]));


nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

figure; set(gcf,'color','w');
subplot(1,2,1);
scatter(d_beh,d_pf,50,'markerfacecolor','k','markeredgecolor','w'); 
xlabel('\DeltaSwitches'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_beh(nan_idx)',d_pf(nan_idx))
text(-38,-0.7,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)


subplot(1,2,2);
scatter(d_meg,d_pf,50,'markerfacecolor','k','markeredgecolor','w'); 
xlabel('\DeltaFC'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_meg(nan_idx),d_pf(nan_idx))
text(-0.018,-0.7,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)


%% correlate baseline
ipharm = 2;
iblock = 1:2;

ifoi = 6;

d_beh = nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1);
d_pf = nanmean(nanmean(b(:,ipharm,iblock),3),2);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7),4));

nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

figure; set(gcf,'color','w');
subplot(1,2,1);
scatter(d_beh,d_pf); 
xlabel('\DeltaSwitches'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_beh(nan_idx)',d_pf(nan_idx))
text(2,8.8,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)

subplot(1,2,2);
scatter(d_meg,d_pf); 
xlabel('\DeltaFC'); ylabel('\DeltaPeak frequency [Hz]');
tp_editplots; axis square; lsline
[r,p]=corr(d_meg(nan_idx),d_pf(nan_idx))
text(0.001,8.8,sprintf('r=%.3f\np=%.3f',r,p),'fontsize',7)

%% topos

ipharm = 2;
iblock = 1:2;

ifoi = 6;

d_beh = nanmean(nanmean(behav_cnt(ipharm,:,iblock),3),1);
d_beh = repmat(d_beh,[274,1]);
d_pf = nanmean(nanmean(pf_all(:,:,ipharm,iblock),4),3);
d_meg = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(fc(:,:,:,ipharm,2,ifoi,iblock),1),2),6),7),4));
d_meg = repmat(d_meg',[274 1]);

nan_idx = ~isnan(d_pf(:)) & ~isnan(d_beh(:));

dd_pf = 100*(nanmean(nanmean(pf_all(:,:,ipharm,iblock),4),3)-nanmean(nanmean(pf_all(:,:,1,iblock),4),3))./nanmean(nanmean(pf_all(:,:,1,iblock),4),3);


% load /home/tpfeffer/pconn/proc/src/pconn_sa_s34_m1_b1_v1.mat
[r,p]=corr(d_beh',d_pf');

pars = [];
figure; set(gcf,'color','white');
pars.scale=[-5 5];
pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.cmap = cmap;
pars.resolution = 300;
showfield_colormap(mean(dd_pf,2),sa.locs_2D,pars);