p_doubledissociation_corr% 
% Decoding using gaussian classifier:
%
% Adrian 05-07-2017
%--------------------------------------------------------------------------

% clear all

ns = 17; % number of subjects
N = 160; % number of nodes

COVar = zeros(N,N,ns); % covariance Atomoxetine-Rest
COVdr = zeros(N,N,ns); % covariance Donepezil-Rest
COVat = zeros(N,N,ns); % covariance Atomoxetine-Task
COVdt = zeros(N,N,ns); % covariance Donepezil-Task


% Get covariances
%-------------------------------------------------------------------------

for subj = 1:ns

% Atomoxetine
%-------------

% Atomoxetine-REST:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series; % X : time-series T-by-N; log(x^2) as in Hipp et al.
x = zscore(X);
COVar(:,:,subj) = cov(x);

% Atomoxetine-TASK:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series; % X : time-series T-by-N; log(x^2) as in Hipp et al.
x = zscore(X);
COVat(:,:,subj) = cov(x);
    
% Donepezil-REST:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series; % X : time-series T-by-N; log(x^2) as in Hipp et al.
T = size(X,1); % number of time steps
x = zscore(X);
COVdr(:,:,subj) = cov(x);

% Donepezil-TASK:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series; % X : time-series T-by-N; log(x^2) as in Hipp et al.
T = size(X,1); % number of time steps
x = zscore(X);
COVdt(:,:,subj) = cov(x);


end



display('decoding...')
% Decode
%-------------------------------------------------------------------------

Class_ar = zeros(1,ns);
Class_dr = zeros(1,ns);
Class_at = zeros(1,ns);
Class_dt = zeros(1,ns);

MU = zeros(1,N);

for subj = 1:ns

trainset = setxor(1:ns,subj); % all subjects except subj.    
COVar_train = mean(COVar(:,:,trainset),3);    
COVdr_train  = mean(COVdr(:,:,trainset),3);    
COVat_train = mean(COVar(:,:,trainset),3);    
COVdt_train  = mean(COVdr(:,:,trainset),3);    
    
% Load single subject data:

% Atomoxetine-REST:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series; % X : time-series T-by-N; log(x^2) as in Hipp et al.
T = size(X,1); % number of time steps
x = zscore(X);

    % decode:
    %-----------------------------------------------------------------
    % Calculate likelihoods and classify:
    groundtruth = 1;  
    hit = classify_data(x,MU,groundtruth,COVar_train,COVdr_train,COVat_train,COVdt_train);
    Class_ar(subj) = hit;

% Atomoxetine-TASK:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series;
x = zscore(X);

    % decode:
    %-----------------------------------------------------------------
    % Calculate likelihoods and classify:
    groundtruth = 3;  
    hit = classify_data(x,MU,groundtruth,COVar_train,COVdr_train,COVat_train,COVdt_train);
    Class_at(subj) = hit;


% Donepezil-REST:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series;
x = zscore(X);

    % decode:
    %-----------------------------------------------------------------
    % Calculate likelihoods and classify:
    groundtruth = 2;  
    hit = classify_data(x,MU,groundtruth,COVar_train,COVdr_train,COVat_train,COVdt_train);
    Class_dr(subj) = hit;
    
% Donepezil-TASK:
file = ; % file containing the time series of subject subj;
load(file);
X = time_series;
x = zscore(X);

    % decode:
    %-----------------------------------------------------------------
    % Calculate likelihoods and classify:
    groundtruth = 4;  
    hit = classify_data(x,MU,groundtruth,COVar_train,COVdr_train,COVat_train,COVdt_train);
    Class_dt(subj) = hit;

end

% percentage of correct classifications:
PCar = 100*sum(Class_a)/length(Class_a); 
PCdr = 100*sum(Class_d)/length(Class_d); 
PCat = 100*sum(Class_a)/length(Class_a); 
PCdt = 100*sum(Class_d)/length(Class_d); 
    

% get chance level(confidence intervals)
% ( = expected performance of a random classifier)
%--------------------------------------------------------------------------
 p=1/4;
 signif=0.05;
 P=zeros(1,ns);
 for k=1:ns
     P(k)=nchoosek(ns,k)*p^k*(1-p)^(ns-k); %binomial distribution
 end
 for k=1:ns
     S=sum(P(k:ns));
     if S<0.05
         break
     end        
 end
Chance=k-1;
Perc95=Chance/ns*100;
 for k=1:ns
     S=sum(P(k:ns));
     if S<0.95
         break
     end   
 end
Chance=k-1;
Perc05=Chance/ns*100;
        
        
figure
plot([.2 4.8],[1 1]*p*100,'k','linewidth',2)
hold on
plot([.2 4.8],[1 1]*Perc95,'k--')
plot([.2 4.8],[1 1]*Perc05,'k--')
bar(1:4,[PCar,PCat,PCdr,PCdt],'barwidth',.9)
xtick = {'Atom. Rest';'Atom. Task';'Done. Rest';'Done. Task'};
set(gca,'xtick',1:4,'xticklabel',xtick)
ylabel('Correct classification (%)')

