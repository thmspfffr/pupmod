%% STATISTICS
% pupmod_all_src_powcorr_map_permtest

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

clear permdat_cnt1
clear permdat_cnt2
clear permdat_res1
clear permdat_res2
v = 12;
outdir   = '/home/tpfeffer/pupmod/proc/conn/';

addpath ~/pconn/matlab/
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

f = fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
aal_labels = textscan(f,'%d %s %d','headerlines',0);


load aalmask_grid_coarse.mat

if v == 1
  for i = 1 : size(cleandat,1)
    for j = 1 : size(cleandat,1)

      idx = find(aalgrid.mask==i,1,'first');
      lab = aalgrid.labels(idx);
      idx = find(strcmp(aal_labels{2},lab));

      jdx = find(aalgrid.mask==j,1,'first');
      lab = aalgrid.labels(jdx);
      jdx = find(strcmp(aal_labels{2},lab));


      p(i,j,:,:,:,:) = cleandat(idx,jdx,:,:,:,:);

    end
  end

  cleandat = p;
  
  for ifoi = 1:13
    %     ifoi

    s_fc(:,:,:,:,1,ifoi) = cleandat(:,:,:,:,1,ifoi);
    s_fc(:,:,:,:,2,ifoi) = cleandat(:,:,:,:,2,ifoi);

  end
else
  for ifoi = 6:7
    %     ifoi

    s_fc(:,:,:,:,1,ifoi) = cleandat(:,:,:,:,1,ifoi);
    s_fc(:,:,:,:,2,ifoi) = cleandat(:,:,:,:,2,ifoi);

  end
end

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')

nperm = 10000; alp = 0.05;

par.subs = 100;
par.allperms = nperm/par.subs;

if ~exist(sprintf('~/pupmod/proc/pupmod_src_powcorr_map_permtest_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  all_idx2 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_map_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v),'all_idx1','all_idx2');
else
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_permtest_map_perms_subs%d_nperm%d_v%d.mat',par.subs,nperm,v));
end

dat_cnt1 = s_fc(:,:,:,[1 2],2,:);
dat_res1 = s_fc(:,:,:,[1 2],1,:);
dat_cnt2 = s_fc(:,:,:,[1 3],2,:); clear s_cnt
dat_res2 = s_fc(:,:,:,[1 3],1,:); clear s_res

taskvsrest(:,:,:,1,:) = dat_res1(:,:,:,1,:);
taskvsrest(:,:,:,2,:) = dat_cnt1(:,:,:,1,:);

for iperm = 1 : par.allperms
  
  if ~exist(sprintf([outdir 'pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' outdir sprintf('pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    
    kperm
    
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
      permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);
      
      permdat_res1(:,:,i,1,:) = dat_res1(:,:,i,idx1(i),:);
      permdat_res1(:,:,i,2,:) = dat_res1(:,:,i,idx2(i),:);
      
    end
%     
%     idx1 = all_idx2(:,iiperm);
%     idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      permdat_cnt2(:,:,i,1,:) = dat_cnt2(:,:,i,idx1(i),:);
      permdat_cnt2(:,:,i,2,:) = dat_cnt2(:,:,i,idx2(i),:);
      
      permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
      permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:);
      
    end
    
    % TASK VS REST ------------------------------------------------------
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      
      taskvsrest_perm(:,:,i,1,:) = taskvsrest(:,:,i,idx1(i),:);
      taskvsrest_perm(:,:,i,2,:) = taskvsrest(:,:,i,idx2(i),:);
      
    end
    
    for ifoi = 6 : 7
      for i = 1 :size(s_fc,1)
        
        [h,~,~,s]=ttest(permdat_res1(i,:,:,2,ifoi),permdat_res1(i,:,:,1,ifoi),'dim',3);
        par.n_p_rest_atx_perm(i,ifoi,kperm) = nansum((h.*sign(s.tstat))>0)./(size(permdat_res1,1)-1);
        par.n_n_rest_atx_perm(i,ifoi,kperm) = nansum((h.*sign(s.tstat))<0)./(size(permdat_res1,1)-1);
        
        [h,~,~,s]=ttest(permdat_cnt1(i,:,:,2,ifoi),permdat_cnt1(i,:,:,1,ifoi),'dim',3);
        par.n_p_task_atx_perm(i,ifoi,kperm) = nansum((h.*sign(s.tstat))>0)./(size(permdat_res1,1)-1);
        par.n_n_task_atx_perm(i,ifoi,kperm) = nansum((h.*sign(s.tstat))<0)./(size(permdat_res1,1)-1);
        
        [h,~,~,s]=ttest(permdat_res2(i,:,:,2,ifoi),permdat_res2(i,:,:,1,ifoi),'dim',3);
        par.n_p_rest_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)>0)./(size(permdat_res1,1)-1);
        par.n_n_rest_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)<0)./(size(permdat_res1,1)-1);
        
        
        [h,~,~,s]=ttest(permdat_cnt2(i,:,:,2,ifoi),permdat_cnt2(i,:,:,1,ifoi),'dim',3);
        par.n_p_task_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)>0)./(size(permdat_res1,1)-1);
        par.n_n_task_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)<0)./(size(permdat_res1,1)-1);
        
        [h,~,~,s]=ttest(permdat_res1(i,:,:,2,ifoi)-permdat_res1(i,:,:,1,ifoi),permdat_cnt1(i,:,:,2,2)-permdat_cnt1(i,:,:,1,2),'dim',3);
        par.n_p_context_atx_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)>0)./(size(permdat_res1,1)-1);
        par.n_n_context_atx_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)<0)./(size(permdat_res1,1)-1);
        
        [h,~,~,s]=ttest(permdat_res2(i,:,:,2,ifoi)-permdat_res2(i,:,:,1,ifoi),permdat_cnt2(i,:,:,2,2)-permdat_cnt2(i,:,:,1,2),'dim',3);
        par.n_p_context_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)>0)./(size(permdat_res1,1)-1);
        par.n_n_context_dpz_perm(i,ifoi,kperm) = nansum(h.*sign(s.tstat)<0)./(size(permdat_res1,1)-1);
        
        % TASK VS REST ------------------------------------------------------
        [t_tvsr,~,~,s] = ttest(taskvsrest_perm(i,:,:,2,ifoi),taskvsrest_perm(i,:,:,1,ifoi),'dim',3,'alpha',alp);
        t_tvsr = t_tvsr.*sign(s.tstat); clear s
        
        par.tperm_taskvsrest_p(i,ifoi,kperm) = nansum(t_tvsr>0)./(size(permdat_res1,1)-1);
        par.tperm_taskvsrest_n(i,ifoi,kperm) = nansum(t_tvsr<0)./(size(permdat_res1,1)-1);
        
        %   NUMBER OF ALTERED CONNECTIONS
        % --------------------------
        %         h=ttest(s_fc(:,:,:,2,1),s_fc(:,:,:,1,1),'dim',3);
        %         atx(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
        %
        %         h=ttest(s_fc(:,:,:,2,2),s_fc(:,:,:,1,2),'dim',3);
        %         atx(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
        %
        %         h=ttest(s_fc(:,:,:,3,1),s_fc(:,:,:,1,1),'dim',3);
        %         dpz(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
        %
        %         h=ttest(s_fc(:,:,:,3,2),s_fc(:,:,:,1,2),'dim',3);
        %         dpz(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
        %
        %         context_allconn_emp_atx(ifoi) = atx(ifoi,2)-atx(ifoi,1);
        %         context_allconn_emp_dpz(ifoi) = dpz(ifoi,2)-dpz(ifoi,1);
        %
        % DOUBLE DISSOCIATION ACROSS FREUQUENCIES -----------------------------
        %         doubledissociation_emp = context_allconn_emp_atx-context_allconn_emp_dpz;
      end    
    end
  end
  
  save(sprintf('~/pupmod/proc/pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
  
  try
    pause(randi(3))
    load(sprintf('~/pupmod/proc/pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  catch me
    save(sprintf('~/pupmod/proc/pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v),'par')
  end
  
end


error('!')
%% PLOT STUFF

nperm         = 10000; 
alp           = 0.05;
par.subs      = 100;
par.allperms  = nperm/par.subs;
v             = 12;

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

% EMPIRICAL STUFF
% --------------------------
for ifoi = [6 7]
  for i = 1 : size(cleandat,1)
    i
    [h,~,~,s]=ttest(cleandat(i,:,:,2,1,ifoi),cleandat(i,:,:,1,1,ifoi),'dim',3);
    n_p_rest_atx(i,ifoi) = nansum((h.*sign(s.tstat))>0)./(size(cleandat,1)-1);
    n_n_rest_atx(i,ifoi) = nansum((h.*sign(s.tstat))<0)./(size(cleandat,1)-1);

    [h,~,~,s]=ttest(cleandat(i,:,:,2,2,ifoi),cleandat(i,:,:,1,2,ifoi),'dim',3);
    n_p_task_atx(i,ifoi) = nansum((h.*sign(s.tstat))>0)./(size(cleandat,1)-1);
    n_n_task_atx(i,ifoi) = nansum((h.*sign(s.tstat))<0)./(size(cleandat,1)-1);

    [h,~,~,s]=ttest(cleandat(i,:,:,3,1,ifoi),cleandat(i,:,:,1,1,ifoi),'dim',3);
    n_p_rest_dpz(i,ifoi) = nansum((h.*sign(s.tstat))>0)./(size(cleandat,1)-1);
    n_n_rest_dpz(i,ifoi) = nansum((h.*sign(s.tstat))<0)./(size(cleandat,1)-1);

    [h,~,~,s]=ttest(cleandat(i,:,:,3,2,ifoi),cleandat(i,:,:,1,2,ifoi),'dim',3);
    n_p_task_dpz(i,ifoi) = nansum((h.*sign(s.tstat))>0)./(size(cleandat,1)-1);
    n_n_task_dpz(i,ifoi) = nansum((h.*sign(s.tstat))<0)./(size(cleandat,1)-1);

% [h,~,~,s]=ttest(permdat_res1(i,:,:,2,ifoi)-permdat_res1(i,:,:,1,ifoi),permdat_cnt1(i,:,:,2,2)-permdat_cnt1(i,:,:,1,2),'dim',3);
% par.n_p_context_atx(i,ifoi) = nansum(h.*sign(s.tstat)>0)./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
% par.n_n_context_atx(i,ifoi) = nansum(h.*sign(s.tstat)<0)./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
% 
% [h,~,~,s]=ttest(permdat_res2(i,:,:,2,ifoi)-permdat_res2(i,:,:,1,ifoi),permdat_cnt2(i,:,:,2,2)-permdat_cnt2(i,:,:,1,2),'dim',3);
% par.n_p_context_dpz(i,ifoi) = nansum(h.*sign(s.tstat)>0)./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
% par.n_n_context_dpz(i,ifoi) = nansum(h.*sign(s.tstat)<0)./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));

% TASK VS REST ------------------------------------------------------
    [t_tvsr,~,~,s] = ttest(cleandat(i,:,:,1,2,ifoi),cleandat(i,:,:,1,1,ifoi),'dim',3,'alpha',alp);
    t_tvsr = t_tvsr.*sign(s.tstat); clear s

    tperm_taskvsrest_p(i,ifoi) = nansum(t_tvsr>0)./(size(cleandat,1)-1);
    tperm_taskvsrest_n(i,ifoi) = nansum(t_tvsr<0)./(size(cleandat,1)-1);

	
  end
end
%%
for iperm = 1 : par.allperms
     
  load(sprintf('~/pupmod/proc/pupmod_src_powcorr_map_permtest_iperm%d_nperm%d_v%d.mat',iperm,nperm,v))
  
  rest_atx_p_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_p_rest_atx_perm;
  rest_atx_n_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_n_rest_atx_perm;  
  task_atx_p_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_p_task_atx_perm;   
  task_atx_n_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_n_task_atx_perm;   
  
  rest_dpz_p_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_p_rest_dpz_perm;
  rest_dpz_n_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_n_rest_dpz_perm;  
  task_dpz_p_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_p_task_dpz_perm;   
  task_dpz_n_perm(:,:,par.subs*(iperm-1)+1:par.subs*iperm) = par.n_n_task_dpz_perm;   
  
  
end

max_rest_atx_p = squeeze(max(rest_atx_p_perm));
max_rest_atx_n = squeeze(max(rest_atx_n_perm));
max_task_atx_p = squeeze(max(task_atx_p_perm));
max_task_atx_n = squeeze(max(task_atx_n_perm));

max_rest_dpz_p = squeeze(max(rest_dpz_p_perm));
max_rest_dpz_n = squeeze(max(rest_dpz_n_perm));
max_task_dpz_p = squeeze(max(task_dpz_p_perm));
max_task_dpz_n = squeeze(max(task_dpz_n_perm));

%% OBTAIN CORRECTED P-VALUES

for ifoi = [6 7]
  for i = 1 : size(cleandat,1)

    p_p_rest_atx(i,ifoi) = 1-sum(n_p_rest_atx(i,ifoi) > max_rest_atx_p(ifoi,:))/nperm;
    p_n_rest_atx(i,ifoi) = 1-sum(n_n_rest_atx(i,ifoi) > max_rest_atx_n(ifoi,:))/nperm;

    p_p_task_atx(i,ifoi) = 1-sum(n_p_task_atx(i,ifoi) > max_task_atx_p(ifoi,:))/nperm;
    p_n_task_atx(i,ifoi) = 1-sum(n_n_task_atx(i,ifoi) > max_task_atx_n(ifoi,:))/nperm;

    p_p_rest_dpz(i,ifoi) = 1-sum(n_p_rest_dpz(i,ifoi) > max_rest_dpz_p(ifoi,:))/nperm;
    p_n_rest_dpz(i,ifoi) = 1-sum(n_n_rest_dpz(i,ifoi) > max_rest_dpz_n(ifoi,:))/nperm;

    p_p_task_dpz(i,ifoi) = 1-sum(n_p_task_dpz(i,ifoi) > max_task_dpz_p(ifoi,:))/nperm;
    p_n_task_dpz(i,ifoi) = 1-sum(n_n_task_dpz(i,ifoi) > max_task_dpz_n(ifoi,:))/nperm;

  end
end

%% IMAGE RESULTS

ifoi = 7;
% var2plot = n_p_task_atx(:,ifoi).*(p_p_task_atx(:,ifoi)<0.025);
var2plot = n_n_rest_dpz(:,ifoi).*(p_n_rest_dpz(:,ifoi)<0.05);


if ~exist('sa_meg_template','var')
  load sa_meg_template;
end

cmap = cbrewer('seq', 'Blues', 256,'pchip'); %cmap(end:-1:1,:);
% cmap = inferno;
% cmap(1:100,:)=0.98*ones(100,3);

mri   = sa_meg_template.mri;

vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = .75;

z2 = spatfiltergauss(var2plot(:),g1,dd,g2);

para = [] ;
para.colorlimits = [min(var2plot) max(var2plot)];

viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];

para.filename = sprintf('~/pupmod/plots/pupmod_all_src_powcorr_map_permtest_f%d_v%d.png',ifoi,v)

tp_showsource(z2,cmap,sa_meg_template,para);

%%








