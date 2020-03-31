%%

load ~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v15.mat
load ~/pmod/matlab/gene_values.mat
idx = tp_aalgenemask(90);

a1 = mean(gdat(idx,find(cell2mat(lab(:,2))==1)),2);
a2 = mean(gdat(idx,find(cell2mat(lab(:,2))==2)),2);
b = mean(gdat(idx,find(cell2mat(lab(:,2))==3)),2);
m = mean(gdat(idx,find(cell2mat(lab(:,2))==7)),2);
n = mean(gdat(idx,find(cell2mat(lab(:,2))==6)),2);
d1 = mean(gdat(idx,find(cell2mat(lab(:,2))==4)),2);
d2 = mean(gdat(idx,find(cell2mat(lab(:,2))==5)),2);

%%
dpz = squeeze(nanmean(nanmean(cleandat(:,:,:,3,1,7),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,:,1,1,7),1),2));
atx = squeeze(nanmean(nanmean(cleandat(:,:,:,2,2,6),1),2))-squeeze(nanmean(nanmean(cleandat(:,:,:,1,2,6),1),2));

squeeze(nanmean(cleandat(:,:,isubj,2,2,6)-cleandat(:,:,isubj,1,2,6),2))


for isubj = 1 : 28
  
  atx = squeeze(nanmean(cleandat(:,:,isubj,2,2,6)-cleandat(:,:,isubj,1,2,6),2));

  corr_NAA1(isubj)  = corr(atx,a1);
  corr_NAA2(isubj)  = corr(atx,a2);
  corr_NAB(isubj)   = corr(atx,b);
  
  corr_D1(isubj)    = corr(atx,d1);
  corr_D2(isubj)    = corr(atx,d2);
  corr_AChM(isubj) 	= corr(atx,m);
  corr_AChN(isubj) 	= corr(atx,n);
  
end

[~,p_NAA1]=ttest(corr_NAA1);
[~,p_NAA2]=ttest(corr_NAA2);
[~,p_NAB]=ttest(corr_NAB);
[~,p_D1]=ttest(corr_D1);
[~,p_D2]=ttest(corr_D2);
[~,p_AChM]=ttest(corr_AChM);
[~,p_AChN]=ttest(corr_AChN);

fprintf('------------------------------\n')
fprintf('Correlations Atomoxetine\n')
fprintf('------------------------------\n')

fprintf('NA/A1: p = %.2f\n',p_NAA1)
fprintf('NA/A2: p = %.2f\n',p_NAA2)
fprintf('NA/B: p = %.2f\n',p_NAB)
fprintf('DA/D1/5: p = %.2f\n',p_D1)
fprintf('DA/D2/3/4: p = %.2f\n',p_D2)
fprintf('ACh/M: p = %.2f\n',p_AChM)
fprintf('ACh/N: p = %.2f\n',p_AChN)


for isubj = 1 : 28
  
  dpz = squeeze(nanmean(cleandat(:,:,isubj,3,1,7)-cleandat(:,:,isubj,1,1,7),2));

  corr_NAA1(isubj)  = corr(dpz,a1);
  corr_NAA2(isubj)   = corr(dpz,a2);
  corr_NAB(isubj)   = corr(dpz,b);
  
  corr_D1(isubj)    = corr(dpz,d1);
  corr_D2(isubj)    = corr(dpz,d2);
  corr_AChM(isubj) 	= corr(dpz,m);
  corr_AChN(isubj) 	= corr(dpz,n);
  
end

[~,p_NAA1]=ttest(corr_NAA1);
[~,p_NAA2]=ttest(corr_NAA2);
[~,p_NAB]=ttest(corr_NAB);
[~,p_D1]=ttest(corr_D1);
[~,p_D2]=ttest(corr_D2);
[~,p_AChM]=ttest(corr_AChM);
[~,p_AChN]=ttest(corr_AChN);

fprintf('------------------------------\n')
fprintf('Correlations Donepezil\n')
fprintf('------------------------------\n')

fprintf('NA/A1: p = %.2f\n',p_NAA1)
fprintf('NA/A2: p = %.2f\n',p_NAA2)
fprintf('NA/B: p = %.2f\n',p_NAB)
fprintf('DA/D1/5: p = %.2f\n',p_D1)
fprintf('DA/D2/3/4: p = %.2f\n',p_D2)
fprintf('ACh/M: p = %.2f\n',p_AChM)
fprintf('ACh/N: p = %.2f\n',p_AChN)




