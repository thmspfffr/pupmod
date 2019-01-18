% fooling around with wavelets

% calculates cross-spectrum and coherence based on a wavelet using
% a Hanning window for given frequeny
%
% usage:  
% [cs coh]=data2cs_wavelet(data,segleng,segshift,epleng,f,fsample);
% input: 
% data:   NxM matrix for N time points and M channels
% segleng:  length of segment (in samples)
% segshift:  shift of segments (in samples)
% epleng: length of epoch (or trial) (in samples)
% f:   frequency of interest (in Hz)
% fsample: sampling frequeny (in Hz)
%
% outpot: 
% cs: cross-spectrum 
% coh: coherency (complex)
% ss: the complex wavelet


segleng = 400;
segshift = 400;
fsample = 400;

f = 10;


%%

load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s4_m1_b1_v6.mat'));
data_low.trial{1} =  data_low.trial{1} + data_hi.trial{1};
[mydata,epleng] = megdata2mydata(data_low); clear data_low data_hi

data = mydata(:,1);
%%

nn=(1:segleng)'-segleng/2;
mywin=hanning(segleng);
s1=cos(nn*f*2*pi/fsample).*mywin;
s2=sin(nn*f*2*pi/fsample).*mywin;
ss=s1-sqrt(-1)*s2;
ssout=ss;

[n nchan]=size(data);
ss=repmat(ss,1,nchan);

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

cs=zeros(nchan,nchan);
coh=cs;


kk=0;
for i=1:nep;
   dloc=data((i-1)*epleng+1:i*epleng,:);
    for j=1:nseg
        kk=kk+1;
        dloc2=dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        dataf=transpose(sum(dloc2.*ss));
        
    end
end




