
%% CREATE CIRCULAR PLOT
% -------------------------------
% In AAL space
addpath(genpath('~/Downloads/3d_plots/'))
v = 1;
cleandat = pupmod_loadpowcorr(v,1);

para = [];
para.N = 90;
para.transfer = 'to_bcn';

for isubj = 1 : 28
  isubj
  for im = 1 : 3
    for icontext = 1 : 2
      for ifoi = 1 : 13
        
        fc(:,:,isubj,im,icontext,ifoi) = tp_match_aal(para,cleandat(:,:,isubj,im,icontext,ifoi));
        
      end
    end
  end
end

cleandat = fc; clear fc;

mask = logical(tril(ones(76,76),-1));
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

nperm = 2000;

cleandat = cleandat(include_bcn,include_bcn,:,:,:,:);

for iperm = 1 : nperm
    
    % within subjects permutation test 
    fprintf('Perm #%d\n',iperm);
    
    idx1 = randi(2,[28 1]);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(:,:,i,1) = cleandat(:,:,i,idx1(i),2,6);
      permdat(:,:,i,2) = cleandat(:,:,i,idx2(i),2,6);     
    end
    
    [~,~,~,s]=ttest(permdat(:,:,:,2),permdat(:,:,:,1),'dim',3);
    maxt_atx(iperm) = max(s.tstat(mask));
        
    idx1 = randi(2,[28 1]); 
    idx2 = 3-idx1; 
    idx1(idx1==2)=3;
    idx2(idx2==2)=3;
    
    for i = 1 : length(idx1)
      permdat(:,:,i,1) = cleandat(:,:,i,idx1(i),1,7);
      permdat(:,:,i,2) = cleandat(:,:,i,idx2(i),1,7);     
    end
    
    [~,~,~,s]=ttest(permdat(:,:,:,2),permdat(:,:,:,1),'dim',3);
    maxt_dpz(iperm) = min(s.tstat(mask));
  
end
%%
% STATS
[h,p,~,s] = ttest(cleandat(:,:,:,2,2,6),cleandat(:,:,:,1,2,6),'dim',3);
mean_change = nanmean(cleandat(:,:,:,2,2,6)-cleandat(:,:,:,1,2,6),3);
conn = s.tstat>prctile(maxt_atx,95);

brain3d_sphere(s.tstat,[prctile(maxt_atx,97.5)],[0 3],[1 5],1,include_bcn)

[h,p,~,s] = ttest(cleandat(:,:,:,3,1,7),cleandat(:,:,:,1,1,7),'dim',3);
mean_change = nanmean(cleandat(:,:,:,3,1,7)-cleandat(:,:,:,1,1,7),3);

brain3d_sphere(abs(s.tstat),[abs(prctile(maxt_dpz,5))],[0 3],[1 5],1,include_bcn)


% STATS
% [h,p,~,s] = ttest(cleandat(include_bcn,include_bcn,:,2,2,6),cleandat(include_bcn,include_bcn,:,1,2,6),'dim',3);
% mean_change = nanmean(cleandat(include_bcn,include_bcn,:,2,2,6)-cleandat(include_bcn,include_bcn,:,1,2,6),3);
% % mean_change = (mean_change-mean(mean_change(mask)))./std(mean_change(mask))
% conn = s.tstat>prctile(maxt_atx,95);
% 
% mean_change = 2 .* ((mean_change - min(mean_change(mask))) ./ (max(mean_change(mask))-min(mean_change(mask))))-1;
% mean_change = mean_change.*conn,labels;
% mean_change(mean_change==0)=nan;
% labels = aal_labels_bcn; 
% labels = labels(include_bcn);
% labels(sum(conn)==0) = {''};
% 
% % figure; set(gcf,'color','w')
% h=schemaball(mean_change,labels,[1 1 1],[0 0 0]);
% 
% print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_schemaball_atx.pdf'));
% %%
% 
% [h,p,~,s] = ttest(cleandat(include_bcn,include_bcn,:,3,1,7),cleandat(include_bcn,include_bcn,:,1,1,7),'dim',3);
% conn = s.tstat<prctile(maxt_dpz,5);
% mean_change = nanmean(cleandat(include_bcn,include_bcn,:,3,1,7)-cleandat(include_bcn,include_bcn,:,1,1,7),3);
% mean_change = 2 .* ((mean_change - min(mean_change(mask))) ./ (max(mean_change(mask))-min(mean_change(mask))))-1;
% mean_change = mean_change.*conn,labels;
% mean_change(mean_change==0)=nan;
% 
% labels = aal_labels_bcn; 
% labels = labels(include_bcn);
% labels(sum(conn)==0) = {''};
% 
% schemaball(abs(mean_change),labels,[1 1 1]);
% 
% print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_all_src_schemaball_dpz.pdf'));

