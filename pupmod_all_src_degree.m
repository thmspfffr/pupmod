%% pupmod_task_src_powcorr

clear

v = 12;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

outdir = '~/pupmod/proc/conn/';
  
ord = pconn_randomization;

for ifoi = 6:7
  
  for isubj = SUBJLIST
    disp(isubj)
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      for iblock = 1 : 2
        clear tmp
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
% 
        s_fc(:,:,isubj,m,1,ifoi,iblock) = powcorr;
        
        load(sprintf([outdir 'pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,v));
%         
        s_fc(:,:,isubj,m,2,ifoi,iblock) = powcorr;

      end
    end
  end
end

s_fc = nanmean(s_fc(:,:,SUBJLIST,:,:,:,:),7);
  
error('!')

%%
for icond = 1 : 1
  for ifoi = 1:13
    
    j = 1:size(s_fc,1);


    for il = 1 : size(s_fc,1)

      jj = j(j~=il);

      for jl = 1 : size(s_fc,1)

        jjj = j(j~=jl);

        fc_tmp = s_fc(:,:,:,1:3,icond,ifoi);

        fprintf('Computing location %d %d ...\n',il,jl);

        x = squeeze(fc_tmp(il,jl,:,:));

        x_ref1 = squeeze(nanmean(fc_tmp(il,jj,:,:),4));
        x_ref2 = squeeze(nanmean(fc_tmp(jjj,jl,:,:),4));

        m1 = mean(x_ref1); s1 = std(x_ref1);
        m2 = mean(x_ref2); s2 = std(x_ref2);

        z11 = (x(:,1)' - m1)./s1;
        z12 = (x(:,1)' - m2)./s2;

        z21 = (x(:,2)' - m1)./s1;
        z22 = (x(:,2)' - m2)./s2;

        z31 = (x(:,3)' - m1)./s1;
        z32 = (x(:,3)' - m2)./s2;

        [~,p11] = ttest(z11);
        [~,p12] = ttest(z12);

        [~,p21] = ttest(z21);
        [~,p22] = ttest(z22);

        [~,p31] = ttest(z31);
        [~,p32] = ttest(z32);

        th11(jl,:) = p11 < 0.05/2;
        th12(jl,:) = p12 < 0.05/2;
        th21(jl,:) = p21 < 0.05/2;
        th22(jl,:) = p22 < 0.05/2;
        th31(jl,:) = p31 < 0.05/2;
        th32(jl,:) = p32 < 0.05/2;

      end

      th11 = th11(jj,:);
      th12 = th12(jj,:);
      th21 = th21(jj,:);
      th22 = th22(jj,:);
      th31 = th31(jj,:);
      th32 = th32(jj,:);

      % any connection significant?
      th1(il,jj,:) = (th11 + th12) > 0; clear th11 th12
      th2(il,jj,:) = (th21 + th22) > 0; clear th21 th22
      th3(il,jj,:) = (th31 + th32) > 0; clear th31 th32

    end

    k(1,ifoi,icond) = sum(sum(th1))/8010
    k(2,ifoi,icond) = sum(sum(th2))/8010
    k(3,ifoi,icond) = sum(sum(th3))/8010

  end
end

%%

figure; 
subplot(1,2,1); hold on

plot(k(1,:,1)'.*100,'color',[.7 .7 .7],'linewidth',2)
plot(k(2,:,1)'.*100,'color',[1 0.5 0.2],'linewidth',2)
plot(k(3,:,1)'.*100,'color',[0.2 0.5 1],'linewidth',2)

xlabel('Frequency [Hz]'); ylabel('Degree [%]')

set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

box off; axis square; axis([1 14 0 35])

subplot(1,2,2); hold on

plot(k(1,:,2)'.*100,'color',[.7 .7 .7],'linewidth',2)
plot(k(2,:,2)'.*100,'color',[1 0.5 0.2],'linewidth',2)
plot(k(3,:,2)'.*100,'color',[0.2 0.5 1],'linewidth',2)

xlabel('Frequency [Hz]'); ylabel('Degree [%]')

set(gca,'tickdir','out','xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])

box off; axis square; axis([1 14 0 35])

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_degree.pdf'))

%% PERMUTATION TEST


% 
% for icond = 1 : 2
% for ifoi = 1:13
% j = 1:90;
% 
% 
% for il = 1 : size(s_fc,1)
%   
%   jj = j(j~=il);
%   
%   for jl = 1 : size(s_fc,1)
%     
%     jjj = j(j~=jl);
%     
%     fc_tmp = s_fc(:,:,:,1:3,icond,ifoi);
%     
%     fprintf('Computing location %d %d ...\n',il,jl);
%     
%     x = squeeze(fc_tmp(il,jl,:,:));
%     
%     x_ref1 = squeeze(nanmean(fc_tmp(il,jj,:,:),4));
%     x_ref2 = squeeze(nanmean(fc_tmp(jjj,jl,:,:),4));
%     
%     m1 = mean(x_ref1); s1 = std(x_ref1);
%     m2 = mean(x_ref2); s2 = std(x_ref2);
%     
%     z11 = (x(:,1)' - m1)./s1;
%     z12 = (x(:,1)' - m2)./s2;
%     
%     z21 = (x(:,2)' - m1)./s1;
%     z22 = (x(:,2)' - m2)./s2;
%     
%     z31 = (x(:,3)' - m1)./s1;
%     z32 = (x(:,3)' - m2)./s2;
%        
%     [~,p11] = ttest(z11);
%     [~,p12] = ttest(z12);
%     
%     [~,p21] = ttest(z21);
%     [~,p22] = ttest(z22);
%     
%     [~,p31] = ttest(z31);
%     [~,p32] = ttest(z32);
%     
%     th11(jl,:) = p11 < 0.01/2;
%     th12(jl,:) = p12 < 0.01/2;
%     th21(jl,:) = p21 < 0.01/2;
%     th22(jl,:) = p22 < 0.01/2;
%     th31(jl,:) = p31 < 0.01/2;
%     th32(jl,:) = p32 < 0.01/2;
%     
%   end
%   
%   th11 = th11(jj,:);
%   th12 = th12(jj,:);
%   th21 = th21(jj,:);
%   th22 = th22(jj,:);
%   th31 = th31(jj,:);
%   th32 = th32(jj,:);
%   
%   % any connection significant?
%   th1(il,jj,:) = (th11 + th12) > 0; clear th11 th12
%   th2(il,jj,:) = (th21 + th22) > 0; clear th21 th22
%   th3(il,jj,:) = (th31 + th32) > 0; clear th31 th32
%     
% end
% 
% k(1,ifoi,icond) = sum(sum(th1))/8010
% k(2,ifoi,icond) = sum(sum(th2))/8010
% k(3,ifoi,icond) = sum(sum(th3))/8010
% 
% end
% 
% end

for ifoi = 1 : 13
  tmp = squeeze(s_fc(:,:,:,1,:,ifoi));
  
  [h,~,~,s] = ttest(tmp(:,:,:,2),tmp(:,:,:,1),'dim',3);
  
  n_p(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)>0))/4005;
  n_n(ifoi) = sum(sum(triu((h.*sign(s.tstat)),1)<0))/4005;
  
end



for iperm = 1 : nperm
  
  
  
  

  
  

