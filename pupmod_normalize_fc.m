function norm_fc = pupmod_normalize_fc(fc)


if size(fc,7)==1
for isubj = 1 : size(fc,3)
    
      tmp = squeeze(fc(:,:,isubj,:,:,:));
      conc = cat(3,tmp(:,:,:,1,:),tmp(:,:,:,2,:));
      
      tmp = zscore(conc,0,3);
      
      norm_fc(:,:,isubj,1:3,1,:) = tmp(:,:,1:3,:);
      norm_fc(:,:,isubj,1:3,2,:) = tmp(:,:,4:6,:);
      
     
  
end

else
    norm_fc = zeros(size(fc));
   for isubj = 1 : size(fc,3)
    
      tmp = squeeze(fc(:,:,isubj,:,:,:,:));
      conc = cat(3,tmp(:,:,:,1,:,1),tmp(:,:,:,2,:,1),tmp(:,:,:,1,:,2),tmp(:,:,:,2,:,2));
      
      tmp = (conc-min(conc,[],3))./(max(conc,[],3)-min(conc,[],3));
      
      norm_fc(:,:,isubj,1:3,1,:,1) = tmp(:,:,1:3,:,:);
      norm_fc(:,:,isubj,1:3,2,:,1) = tmp(:,:,4:6,:,:);
      norm_fc(:,:,isubj,1:3,1,:,2) = tmp(:,:,7:9,:,:);
      norm_fc(:,:,isubj,1:3,2,:,2) = tmp(:,:,10:12,:,:);
      clear tmp

   end
end
      