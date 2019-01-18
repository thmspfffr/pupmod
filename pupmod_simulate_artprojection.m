%% SIMULTATE PROJECTION ARTIFACTS
regression = 1;allcond = 1;

NSUBJ = 28;
clc; clear fc cnt1 cnt2 true cleandat
% close all; 
nsim = 100;
par = [1.2 1.5];

for isim = 1 : nsim
  %   clc
  fprintf('isim%d...\n............\n',isim)
  
  cleandat = zeros(90,90,NSUBJ,2);
  
  n1 = rand(NSUBJ,1);
  nuis1 = [n1 ((rand(NSUBJ,1)+0.5)*par(1)).*n1];
  n2 = rand(NSUBJ,1);
  nuis2 = [n2 ((rand(NSUBJ,1)+0.5)*par(2)).*n2];
  
  base1 = rand(90,90,NSUBJ);
  base2 = 0.8*rand(90,90,NSUBJ);
  
  [h,~,~,t]=ttest(base2,base1,'dim',3);
  true(isim,1) = sum(sum((t.tstat>0).*h))/8100;
  true(isim,2) = sum(sum((t.tstat<0).*h))/8100;
  
  fprintf('True: P = %.3f | N = %.3f\n',true(isim,1),true(isim,2))
  
  fc(:,:,:,1) = base1+permute(repmat(nuis1(:,1),[1 90 90]),[2 3 1])+permute(repmat(nuis2(:,1),[1 90 90]),[2 3 1]);
  fc(:,:,:,2) = base2+permute(repmat(nuis1(:,2),[1 90 90]),[2 3 1])+permute(repmat(nuis2(:,2),[1 90 90]),[2 3 1]);
  
  [h,~,~,t]=ttest(fc(:,:,:,2),fc(:,:,:,1),'dim',3);
  cnt1(isim,1) = sum(sum((t.tstat>0).*h))/8100;
  cnt1(isim,2) = sum(sum((t.tstat<0).*h))/8100;
  
  fprintf('Artf: P = %.3f | N = %.3f\n',cnt1(isim,1),cnt1(isim,2))
  
  if allcond
    for i = 1 : 90
      % 	fprintf('Cond %d Session %d freq %d node %d...\n',icont,im,ifoi,i)
      for j = 1 : 90
        siz = size(squeeze(fc(i,j,:,:)));
        
        dat = reshape(squeeze(fc(i,j,:,:)),[prod(siz) 1]);
        
        if ~regression
          nuisance_var = [zscore(nuis1(:)) zscore(nuis2(:))];
          tmp = decorrgs(nuisance_var,size(nuisance_var,2));
          
          refs = zscore(tmp)./norm(zscore(tmp));
          
          for inui = 1 : size(nuisance_var,2)
            dat = (dat - (dat'*refs(:,inui))*refs(:,inui));
          end        
        else
          tmp = [zscore(nuis1(:)) zscore(nuis2(:))];
          [~,~,dat]=regress(dat,[tmp(:,1) tmp(:,2)]);
        end
        cleandat(i,j,:,:) = reshape(dat,siz);
      end
    end
  else
    for im = 1 : 2
      for i = 1 : 90
        for j = 1 : 90
          
          siz = size(squeeze(fc(i,j,:,im)));
          dat = reshape(squeeze(fc(i,j,:,im)),[prod(siz) 1]);
          
          if ~regression
            siz = size(squeeze(fc(i,j,:,im)));
            nuisance_var = [zscore(nuis1(:,im)) zscore(nuis2(:,im))];
            tmp = decorrgs(nuisance_var,size(nuisance_var,2));
            
            refs = zscore(tmp)./norm(zscore(tmp));
            
            for inui = 1 : size(nuisance_var,2)
              dat = (dat - (dat'*refs(:,inui))*refs(:,inui));
            end
          else
            [~,~,dat]=regress(dat,[nuis1(:,im) nuis2(:,im)]);
          end
          cleandat(i,j,:,im) = reshape(dat,siz);
        end
      end
    end
  end
  
  [h,~,~,t]=ttest(cleandat(:,:,:,2),cleandat(:,:,:,1),'dim',3);
  cnt2(isim,1) =  sum(sum((t.tstat>0).*h))/8100;
  cnt2(isim,2) =  sum(sum((t.tstat<0).*h))/8100;
  
  fprintf('Clean: P = %.3f | N = %.3f\n',cnt2(isim,1),cnt2(isim,2))
  
end

%%
figure; set(gcf,'color','w')

m = [mean(true(1:nsim,:)); mean(cnt1(1:nsim,:)); mean(cnt2(1:nsim,:))];
s = [std(true(1:nsim,:),[],1); std(cnt1(1:nsim,:),[],1); std(cnt2(1:nsim,:),[],1)];

% subplot(1,2,1); hold on
bar([0.75 2.75 4.75],m(:,1),0.4,'facecolor','r'); hold on
line([0.75 0.75; 2.75 2.75; 4.75 4.75]',[m(:,1)-s(:,1) m(:,1)+s(:,1)]','color','k')
bar([1.25 3.25 5.25],m(:,2),0.4,'facecolor','b'); hold on
line([1.25 1.25; 3.25 3.25; 5.25 5.25]',[m(:,2)-s(:,2) m(:,2)+s(:,2)]','color','k')

tp_editplots;
set(gca,'XTick',1:2:5,'xTicklabel',{'True', 'Contaminated','Cleaned'})
ylabel('Fraction of altered correlations')

