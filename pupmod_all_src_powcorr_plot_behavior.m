%% pupmod_all_src_powcorr_plot
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

v = 12;

% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
SUBJLIST  = [4 5 6 7 8 9 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

% error('!')

%%

mask = find(triu(ones(400))-eye(400));
fc_all = reshape(cleandat(:,:,:,1:3,1:2,:),[400*400 28 3 2 13]);
fc_all = fc_all(mask,:,:,:,:); %clear cleandat

% load ~/pconn_bttn/proc/behav_counting_upload.mat
load ~/pconn_bttn/proc/behav_pressing_upload.mat

behav = (behav_pressing + behav_counting) / 2;
behav(7,2) = behav_counting(7,2);
behav(27,2) = behav_counting(27,2);
% behav([7,27],:) = []

clear outp


%% CORRELATE CHANGE IN BEHAVIOR WITH CHANGE IN FC
% Preselect only connections that are altered through task
for ifoi = 1:13
  ifoi
  idx1 = ttest(fc_all(:,:,2,2,ifoi),fc_all(:,:,1,2,ifoi),'dim',2,'alpha',0.01);
  idx2 = ttest(fc_all(:,:,2,2,ifoi),fc_all(:,:,2,1,ifoi),'dim',2,'alpha',0.01);

  fc = mean(fc_all(find(idx1),:,2,2,ifoi)-fc_all(find(idx1),:,1,2,ifoi));

  [outp.rr1(:,ifoi) outp.pp1(:,ifoi)] = corr(fc',behav(:,2)-behav(:,1));
  
  fc = mean(fc_all(find(idx1&idx2),:,2,2,ifoi)-fc_all(find(idx1&idx2),:,1,2,ifoi));
  
  [outp.rr2(:,ifoi) outp.pp2(:,ifoi)] = corr(fc',behav(:,2)-behav(:,1));
  
  fc = mean(fc_all(find(idx1),:,2,2,ifoi)-fc_all(find(idx1),:,1,2,ifoi));

  [outp.rr3(:,ifoi) outp.pp3(:,ifoi)] = corr(fc',behav_counting(:,2)-behav_counting(:,1));
  
  fc = mean((fc_all(find(idx1&idx2),:,2,2,ifoi)-fc_all(find(idx1&idx2),:,1,2,ifoi)));
  
  [outp.rr4(:,ifoi) outp.pp4(:,ifoi)] = corr(fc',behav_counting(:,2)-behav_counting(:,1));
  
end

figure; set(gcf,'color','w'); hold on
subplot(1,2,1); hold on
plot(outp.rr1)
plot(outp.rr2)
title('Pressing')
xlabel('Frequencies'); ylabel('Correlation (FC, Behav)'); axis square
axis([0 14 -0.4 0.2])

subplot(1,2,2); hold on
plot(outp.rr3)
plot(outp.rr4)
title('Counting')
xlabel('Frequencies'); ylabel('Correlation (FC, Behav)'); axis square
axis([0 14 -0.4 0.2])


%% CORRELATE CHANGE IN BEHAVIOR WITH CHANGE IN FC
% First find connections that correlate with behavior

for ifoi = 1:13
  ifoi
  for i = 1 : 79800
    i
    [r(i), p(i)] = corr(fc_all(i,:,1,2,ifoi)',behav(:,1));
  end
  idx1 = p<0.05;
  
  for i = 1 : 79800
    i
    [r(i), p(i)] = corr(fc_all(i,:,1,2,ifoi)',behav_counting(:,1));
  end
  fc = mean(fc_all(find(idx1),:,2,2,ifoi)-fc_all(find(idx1),:,1,2,ifoi));

  [outp.rr1(:,ifoi) outp.pp1(:,ifoi)] = corr(fc',behav(:,2)-behav(:,1));
  [outp.rr2(:,ifoi) outp.pp2(:,ifoi)] = corr(fc',behav_counting(:,2)-behav_counting(:,1));

end

figure; set(gcf,'color','w'); hold on
subplot(1,2,1); hold on
plot(outp.rr1)
plot(outp.rr2)
title('Pressing')
xlabel('Frequencies'); ylabel('Correlation (FC, Behav)'); axis square
axis([0 14 -0.4 0.2])


%% load('~/pupmod/proc/pupmod_src_behavcorr_tsk_cnt.mat')

ifoi = 6;
context = 2;
drug = 2; 
alpha = 0.05;

[idx] = find(outp.pp(:,ifoi) < alpha);

fc = reshape(cleandat(:,:,:,drug,context,ifoi),[400*400 28])-reshape(cleandat(:,:,:,1,context,ifoi),[400*400 28]);
fc = fc(mask,:);

[r,p]=corr(mean(fc(idx,:))',behav(:,drug)-behav(:,1))

figure; set(gcf,'color','w')
scatter(mean(fc(idx,:)),behav(:,drug)-behav(:,1))
lsline

axis square
xlabel('Difference in FC (Atx - Pbo)')
ylabel('Difference in behavior (Atx - Pbo)')
tp_editplots


% print(gcf,'dpdf','~/Dropbox/projects/phd/pupmod/plots/plot.pdf')
%% PERMUTAITON TEST
clear r_perm

mask = logical(tril(ones(400,400),-1));
mask = find(triu(ones(400))-eye(400));

fc = reshape(cleandat(:,:,:,1:2,2,6),[400*400 28 2]);
fc = fc(mask,:,:);

for iperm = 1 : 1000
    iperm
    idx1 = randi(2,[28,1]);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat(:,i,1) = fc(:,i,idx1(i));
      permdat(:,i,2) = fc(:,i,idx2(i));
      
    end
    
%     [rr, pp] = corr(squeeze(permdat(:,:,1))',behav(:,1));
%     idx2 = find(pp < alpha);
%     
%     if isempty(idx2)
%         r_perm(iperm) = 0;
%         continue
%     end
    
    d_meg = permdat(:,:,2)-permdat(:,:,1);
    r_perm(iperm)=corr(mean(d_meg(idx,:),1)',behav(:,drug)-behav(:,1));

end
    
    
    
p  =1-sum(r<r_perm(1:1000))/1000
    

    








