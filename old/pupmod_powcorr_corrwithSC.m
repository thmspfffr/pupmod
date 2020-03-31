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