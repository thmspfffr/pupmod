%% pupmod_all_regressartifacts
% regress out artifacts (from pupmod_all_powcorr_peripheral)


clear

% vv = 11; % muscles
% vv =12 ; % heart/muscles
% vv =13 ; % blinks/heart/muscles

v = 23; vv = 23;

outdir = '~/pupmod/proc/conn/';
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

load(['/home/tpfeffer/pconn_cnt/proc/' sprintf('pupmod_all_powcorr_peripheral.mat')])

indiv_subj = 0;

%% load data


if indiv_subj == 1
  
  fc = pupmod_loadpowcorr(v,SUBJLIST,0);

  % PERFORM REGRESSION *WITHIN SUBJECTS*, BUT
  % ACROSS BLOCKS, CONTEXTS AND CONDITONS
  % ---------------------
  for isubj = SUBJLIST
    isubj
       
    fn = sprintf('pupmod_all_regressartifacts_withinsubj_s%d_v%d',isubj,vv);
    if tp_parallel(fn,'~/pupmod/proc/conn/',1,0)
      continue
    end
    
  	res_dat = zeros(size(fc,1),size(fc,1),siz(1), siz(2),25,siz(3));

    for ifoi = 1:25
   ifoi
      siz = size(squeeze(fc(1,1,isubj==SUBJLIST,:,:,ifoi,:)));
      
      art1 = reshape(par.blinks(isubj==SUBJLIST,:,:,:),[prod(siz) 1]);
      art2 = reshape(par.heartrate(isubj==SUBJLIST,:,:,:),[prod(siz) 1]);
      art3 = reshape(par.muscles(isubj==SUBJLIST,:,:,:),[prod(siz) 1]);
      
      art1 = (art1-nanmean(art1))/nanstd(art1);
      art2 = (art2-nanmean(art2))/nanstd(art2);
      art3 = (art3-nanmean(art3))/nanstd(art3);
      
      nuisance_var = [art1 art2 art3];
            
      for i = 1 : size(fc,1)
        
        for j = 1 : size(fc,1)
          
          dat = reshape(squeeze(fc(i,j,isubj==SUBJLIST,:,:,ifoi,:)),[prod(siz) 1]);
          
          if ~any(isnan(dat))
            [~,~,dat]=regress(dat,nuisance_var);
          elseif sum(isnan(dat))>0 && sum(isnan(dat))<12
            idx = find(isnan(dat));
            [~,~,dat(~isnan(dat))]=regress(dat(~isnan(dat)),nuisance_var(~isnan(dat),:));
          else
            
          end
          
          res_dat(i,j,:,:,ifoi,:) = reshape(dat,siz);
          
        end
      end
    end
    
    save(sprintf('~/pupmod/proc/conn/%s.mat',fn),'res_dat','-v7.3')

  end
  
else


  % PERFORM REGRESSION ACROSS SUBJECTS, BLOCKS, CONTEXTS AND CONDITONS
  % ---------------------
  
  fc = pupmod_loadpowcorr(v,SUBJLIST,1);

  for ifoi = 1:25
    ifoi
    fn = sprintf('pupmod_all_regressartifacts_acrosssubj_f%d_v%d',ifoi,vv);
    if tp_parallel(fn,'~/pupmod/proc/conn/',1,0)
      continue
    end
       

    s_fc = fc(:,:,:,:,:,ifoi);
    
    siz = size(squeeze(s_fc(1,1,:,:,:)));
    
    art1 = reshape(nanmean(par.blinks,4),[prod(siz) 1]);
    art2 = reshape(nanmean(par.heartrate,4),[prod(siz) 1]);
    art3 = reshape(nanmean(par.muscles,4),[prod(siz) 1]);
    
    art1 = (art1-nanmean(art1))/nanstd(art1);
    art2 = (art2-nanmean(art2))/nanstd(art2);
    art3 = (art3-nanmean(art3))/nanstd(art3);
    
    nuisance_var = [art1 art2 art3];
    
    for i = 1 : size(s_fc,1)
      for j = 1 : size(s_fc,1)
        
        dat = reshape(squeeze(s_fc(i,j,:,:,:)),[prod(siz) 1]);
        
        if ~any(isnan(dat))
          [~,~,dat]=regress(dat,nuisance_var);
        elseif sum(isnan(dat))>0 && sum(isnan(dat))<prod(siz)
          idx = find(isnan(dat));
          [~,~,dat(~isnan(dat))]=regress(dat(~isnan(dat)),nuisance_var(~isnan(dat),:));
        else
          error('!')
        end
        
        res_dat(i,j,:,:,:) = reshape(dat,siz);
        
      end
    end 
    save(sprintf('~/pupmod/proc/conn/%s.mat',fn),'res_dat','-v7.3')

  end

end

error('!')

%%
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% % if isubj==length(SUBJLIST)
%   if vv==1
%     cleandat = zeros(90,90,28,3,2,13,2,'single');
%   elseif vv>=11
%     cleandat = zeros(400,400,28,3,2,25,2,'single');
%   end
%   
%   for isubj = 1:length(SUBJLIST)
%     isubj
%     fn = sprintf('pupmod_all_regressartifacts_s%d_v%d',SUBJLIST(isubj),vv);
%     load(sprintf('~/pupmod/proc/conn/%s.mat',fn))
%     
%     cleandat(:,:,isubj,:,:,:,:) = single(res_dat);
%   end
%   
fc = res_dat; clear res_dat; fc = permute(fc,[1 2 3 4 6 5]);
  save(sprintf(['~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat'],vv),'fc','-v7.3');
% end
%%
mask    = logical(triu(ones(400,400),1));

for ifoi = 1 :25
  ifoi
  [h,~,~,s]=ttest(squeeze(res_dat(:,:,:,2,2,ifoi)),squeeze(res_dat(:,:,:,1,2,ifoi)),'dim',3);
  pos_atx(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_atx(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  [h,~,~,s]=ttest(squeeze(res_dat(:,:,:,3,2,ifoi)),squeeze(res_dat(:,:,:,1,2,ifoi)),'dim',3);
  pos_dpz(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_dpz(ifoi,2) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
%   
  [h,~,~,s]=ttest(squeeze(res_dat(:,:,:,2,1,ifoi)),squeeze(res_dat(:,:,:,1,1,ifoi)),'dim',3);
  pos_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  [h,~,~,s]=ttest(squeeze(res_dat(:,:,:,3,1,ifoi)),squeeze(res_dat(:,:,:,1,1,ifoi)),'dim',3);
  pos_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
  neg_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
  
%   [h,~,~,s]=ttest(squeeze(nanmean(pow12_res(:,:,:,:,2,ifoi),3)),squeeze(nanmean(pow12_res(:,:,:,:,1,ifoi),3)),'dim',3);
%   pos_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
%   neg_atx(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
%   [h,~,~,s]=ttest(squeeze(nanmean(pow12_res(:,:,:,:,3,ifoi),3)),squeeze(nanmean(pow12_res(:,:,:,:,1,ifoi),3)),'dim',3);
%   pos_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)>0))/sum(mask(:));
%   neg_dpz(ifoi,1) = 100*sum((h(mask)>0)&(s.tstat(mask)<0))/sum(mask(:));
end



