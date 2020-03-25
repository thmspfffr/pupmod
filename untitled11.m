%% NOISE STIMULATION




data = rand(1000,17,28,2);

sign_idx = logical(randi(2,[1000 1])-1); 

data(sign_idx,8:12,1:28,2)=data(sign_idx,8:12,1:28,2)+(rand(sum(sign_idx),5,28)-0.475);

[h,p,~,s]=ttest(mean(data(:,:,:,2),1),mean(data(:,:,:,1),1),'dim',3);

% altered_p = sum(h>0 & s.tstat > 0)/100;
% altered_n = sum(h>0 & s.tstat < 0)/100;


[h,p,~,s]=ttest(mean(data(:,:,:,2),2),mean(data(:,:,:,1),2),'dim',3);

fprintf('Identified: %d\n', sum(h))
fprintf('Real: %d\n', sum(sign_idx))


%%




data = rand(1000,17,28,2);

sign_idx = logical(randi(2,[1000 1])-1); 

[h,p,~,s]=ttest(mean(data(:,:,:,2),1),mean(data(:,:,:,1),1),'dim',3);

% altered_p = sum(h>0 & s.tstat > 0)/100;
% altered_n = sum(h>0 & s.tstat < 0)/100;


[h,p,~,s]=ttest(mean(data(:,p<0.05,:,2),2),mean(data(:,p<0.05,:,1),2),'dim',3);

fprintf('Identified: %d\n', sum(h))
% fprintf('Real: %d\n', sum(sign_idx))






