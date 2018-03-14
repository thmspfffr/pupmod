%% pupmod_mri_kuramoto
% Compute MRI Kuramoto order parameter

load ~/M.mat

% bp = hilbert(tp_bpfilt(squeeze(M(1,1,:,1)),[0.02 0.2],1/2.25,4,'butter'))

bp = [0.01 0.22]

for isubj = 1 : 24
  isubj
  for icond = 1 : 4
    bp_mri(:,:,isubj,icond) = squeeze(M(isubj,icond,:,:));
%     bp_mri(:,:,isubj,icond) = hilbert(tp_bpfilt(squeeze(M(isubj,icond,:,:)),[bp],1/2.25,4,'butter'));
    kura(:,isubj,icond) = abs(sum(exp(1i*angle(bp_mri(:,:,isubj,icond))),2))/90;
    
  end
end
t=squeeze(mean(kura))
[~,p]=ttest(t(:,3),t(:,4))
t=squeeze(std(kura))
[~,p]=ttest(t(:,3),t(:,4))
