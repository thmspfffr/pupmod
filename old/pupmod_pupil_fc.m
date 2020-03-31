addpath ~/pupmod/matlab

% load fc
fc = pupmod_loadpowcorr(12,1);

% load pupil
para.time = [3000 7000];
pup = pp_loadpupil(para);

% third dimension: 1 = rest, 2 = count, 3 = button
% second dim: 1 = placebo, 2 = atomoxetine, 3 = donepezil
pup = squeeze(nanmean(pup,3));
%%
icont = 1;
ipharm = 1;
crit_freq = 9;

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);
for ifoi = 1 : 13
  if ipharm == 1
    fc_atx  = fc(:,:,:,ipharm,icont,ifoi);
    pup_atx = pup(:,ipharm,icont);
  else
    fc_atx  = fc(:,:,:,ipharm,icont,ifoi)-fc(:,:,:,1,icont,ifoi);
    pup_atx = pup(:,ipharm,icont)-pup(:,1,icont);
  end
tmp_fc  = reshape(fc_atx,[400*400 28]);
tmp_pup = repmat(pup_atx',[400*400 1]);
    
% ------------
% compute correlation
% ------------
a = bsxfun(@minus, tmp_fc, mean(tmp_fc,2));
b = bsxfun(@minus, tmp_pup, mean(tmp_pup,2));

a1 = a.^2;
b1 = b.^2;
ab = a.*b;

r = sum(ab,2) ./ sqrt ( sum(a1,2) .* sum(b1,2));
t = r .* sqrt((28-2)./(1-r.^2));
p = 2*(1-tcdf(t,28-2));

t = r .* sqrt((28-2)./(1-r.^2));
s = tcdf(t,28-2);
p = 2*min(s,1-s);

larger(ifoi)  = 100*(sum((t>0)&(p<0.05))./size(t,1));
smaller(ifoi) = 100*(sum((t<0)&(p<0.05))./size(t,1));

r = reshape(r,[400 400]);
t = reshape(t,[400 400]);
p = reshape(p,[400 400]);

if ifoi == crit_freq
  r6 = r;
  p6 = p;
end
% ------------
% 
% for i = 1 : 400
%   for j = 1 : 400
%     [rr(i,j) pp(i,j)] = corr(squeeze(fc_atx(i,j,:)),pup_atx);
%   end
% end
% 

end

foi_range       = unique(round(2.^[1:.5:7]));
figure; set(gcf,'color','w'); 
subplot(3,4,[1 2 3 4]);hold on
plot(smaller,'linewidth',2,'color','b');
plot(larger,'linewidth',2,'color','r');
line([crit_freq crit_freq],[0 max([smaller(:); larger(:)])],'color','k','linestyle',':')
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
tp_editplots
ylabel('Fraction of sign. corr. [%]'); xlabel('Carrier frequency [Hz]')

h = p6 < 0.05;
subplot(3,4,[5 6 9 10]); 
imagesc(r6,[-0.5 0.5]); axis off square
subplot(3,4,[7 8 11 12]); 
imagesc(r6.*h,[-0.5 0.5]); axis off square;
xlabel('Region')
colormap(cmap);