clear all

cd ~/Downloads/

load SC_AAL_90.mat

areas = 90;
C = SC/max(max(SC))*0.2;

time_steps = 600;

a = 0.2;
we = 0.8;
sig = 0.01;
omega = 10; 

threshold = 0;
polarity = 1;

dt = 0.0025;

xs = zeros(time_steps/2,areas);
x = 0.1 * randn(areas,1);
y = 0.1 * randn(areas,1);
nn = 0;

Isubdiag = find(tril(ones(areas),-1));


for t=0:dt:500
    t
 sumax = C*x-diag(C*repmat(x',areas,1));
 sumay = C*y-diag(C*repmat(y',areas,1));

 x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
 y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
 
end

for t=0:dt:5000  %32000
 t
% if t > 500
%     a = 0.1;   
% end
 
 sumax = C*x-diag(C*repmat(x',areas,1));
 sumay = C*y-diag(C*repmat(y',areas,1));

 x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
 y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
 

if mod(t,2)==0
    nn=nn+1;
    xs(nn,:)=x';
    ri=sqrt(x.*x+y.*y);
    rimax=max(ri);
    ri=ri/rimax;
    kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
end

end
 

FC_simul = corrcoef(xs(1:nn,:));


% RUN DFA
f_win = [1 100];
c_win = [0.5 150];

siginfo = nbt_Info;
siginfo.converted_sample_frequency =  5;

t = nbt_doDFA(abs(hilbert(xs)), siginfo, [1 100],[0.5 150],0.5,0,0,[]);
% cc = corrcoef(FC_emp(Isubdiag),FC_simul(Isubdiag));
v = var(abs(hilbert(xs)),[],1);
m = mean(abs(hilbert(xs)));


%%
% -0.2
m = 0.0117
v = 3.78e-05
d = 0.56

% 0
m = 0.271
v = 2.08e-04
d = 0.92

% 0.2
m = 0.44;
v = 8.72e-05;
d = 0.56

save('~/pmod/proc/pmod_hopf_par2.mat','t','v','m')


Corr = cc(1,2);

Corr_sim_mean = mean(mean(FC_simul));

data = xs';
[Peak,Amp] = Threshold_fMRI(data, areas, threshold, polarity);

figure;
subplot(2,1,1)
pcolor(Peak)
subplot(2,1,2)
plot(1:time_steps/2+1, xs(:,1:66))
xlim([0 time_steps/2+1])

