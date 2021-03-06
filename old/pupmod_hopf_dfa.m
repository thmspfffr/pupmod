%% pupmod_src_hopf
% pupmod_hopf_dfa

clear all

% --------------------------------------------------------
% VERSION 6
% --------------------------------------------------------
% v          = 6;
% fsample    = 400;
% para.bif   = -0.75:0.05:0.75;
% para.G     = 1.5;
% para.sig   = .01:.01:.05;
% para.freq  = 10;
% para.ts    = 600; % timesteps
% plt = 0;
% --------------------------------------------------------
% VERSION 8 - same as 6, but no coupling
% --------------------------------------------------------
v          = 8;
fsample    = 400;
para.bif   = -0.75:0.05:0.75;
para.G     = 0;
para.sig   = [.01 .05];
para.freq  = 10;
para.ts    = 600; % timesteps
plt = 0;
% --------------------------------------------------------

% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
% v          = 7;
% fsample    = 400;
% para.bif   = -0.07:0.005:0.07;
% para.G     = 1.5;
% para.sig   = .01:.03:.15;
% para.freq  = 10;
% para.ts    = 600; % timesteps
% plt = 0;
% --------------------------------------------------------


outdir   = '/home/tpfeffer/pupmod/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%% RUN HOPF MODEL

for isig = 1 : length(para.sig)
  for ibi = 1 : length(para.bif)
    for iG = 1 : length(para.G)
      
      if ~exist(sprintf(['~/pupmod/proc/pupmod_hopf_dfa_n%d_b%d_g%d_v%d_processing.txt'],isig,ibi,iG,v))
        system(['touch ' '~/pupmod/proc/' sprintf('pupmod_hopf_dfa_n%d_b%d_g%d_v%d_processing.txt',isig,ibi,iG,v)]);
      else
        continue
      end
      
      clear xs kura
            
      load ~/Downloads/SC_AAL_90.mat
      
      areas       = 90;
      C           = SC/max(max(SC))*0.2;
      
      time_steps  = para.ts;
      
      a           = para.bif(ibi);  % bifurcation param.
      we          = para.G(iG);     % coupling
      sig         = para.sig(isig);           % amp. of noise?
      f           = para.freq;             % freq. = omega
      
      dt          = 1/fsample;
      
      amp         = zeros(time_steps/dt,areas);
      ys          = zeros(time_steps/dt,areas);
      xs          = zeros(time_steps/dt,areas);
      
      kura        = zeros(time_steps/dt,1);
      
      x           = 0.1 * randn(areas,1);
      y           = 0.1 * randn(areas,1);
      nn          = 0;
      omega       = 2*pi*f;
      
      Isubdiag = find(tril(ones(areas),-1));
      
      fprintf('Initializing...\n');
      
      xs = zeros(para.ts/(1/fsample)+1,90);
      ys = zeros(para.ts/(1/fsample)+1,90);
      amp = zeros(para.ts/(1/fsample)+1,90);
      kura = zeros(para.ts/(1/fsample)+1,1);
      
      for t=0:dt:100
        
        sumax = C*x-diag(C*repmat(x',areas,1));
        sumay = C*y-diag(C*repmat(y',areas,1));
        
        x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
        y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
        
      end
      
      if plt
        figure; set(gcf,'color','w');
      end
      fprintf('Running model...\n');
      ts_cnt = 0;
      for t=0:dt:time_steps  %32000
        ts_cnt = ts_cnt+1;
        fprintf('Time step %.3f / %d ...\n',t,time_steps);
        
        sumax = C*x-diag(C*repmat(x',areas,1));
        sumay = C*y-diag(C*repmat(y',areas,1));
        
        %       x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
        %       y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
        x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
        y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
        
        %       if mod(t,2)==0
        nn=nn+1;
        xs(nn,:)=x';
        ys(nn,:)=y';
        
        amp(nn,:) = sqrt(x.*x+y.*y);
        
        %       if plt
        %         if nn >  100
        %           subplot(2,1,1)
        %           plot(xs(nn:-1:nn-100),ys(nn:-1:nn-100),'-')
        %           axis square tight
        %           %       axis([-sqrt(para.bif)-sqrt(para.bif)*0.2 sqrt(para.bif)+sqrt(para.bif)*0.2 -sqrt(para.bif)-sqrt(para.bif)*0.2 sqrt(para.bif)+sqrt(para.bif)*0.2]); grid on;
        %           subplot(2,1,2)
        %           plot(101:-1:1,amp(nn:-1:nn-100,1))
        %           axis square tight
        %           drawnow
        %         end
        %       end
        
        ri=amp(nn,:)';
        rimax=max(ri);
        ri=ri/rimax;
        kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
        
      end
      
      siginfo = nbt_Info;
      siginfo.converted_sample_frequency = 400;
      
      tmp = nbt_doDFA(amp, siginfo, [1 100],[0.5 150],0,0,0,[]);
      
      par.dfa(:,ibi) = tmp.MarkerValues;
      par.pow(:,ibi) = mean(amp);
      par.var(:,ibi) = var(amp);
      
      par.fc(:,:,ibi) = corrcoef(amp(1:nn,:));
      
      save(sprintf('~/pupmod/proc/pupmod_hopf_dfa_sig%d_bif%d_v%d.mat',isig,ibi,v),'par'); clear par amp xs ys
      
      
    end
  end
end


error('!')

%%
for isig = 1 : length(para.sig)
  for ibi = 1 : length(para.bif)

    load(sprintf('~/pupmod/proc/pupmod_hopf_dfa_sig%d_bif%d_v%d.mat',isig,ibi,v));

    dfa(:,ibi,isig)  = par.dfa(:,ibi);
    pow(:,ibi,isig)  = par.pow(:,ibi);
    vr(:,ibi,isig)    = par.var(:,ibi);
    fc(:,:,ibi,isig) = par.fc(:,:,ibi);

  end
end
%% COMPUTE ALL MEASURES ACROSS BIF PARAMETER

isig = 1;

figure_white; hold on; title(sprintf('Noise: %d',isig))

plot(para.bif,zscore(mean(dfa(:,:,isig))),'linewidth',3);
plot(para.bif,zscore(mean(pow(:,:,isig))),'linewidth',3);
plot(para.bif,zscore(mean(vr(:,:,isig))),'linewidth',3);
plot(para.bif,zscore(squeeze(mean(mean(fc(:,:,:,isig),1),2))),'linewidth',3)

xlabel('Bifurcation parameter'); ylabel('DFA/Pow/Var/FC');

legend({'DFA';'Amp';'Var';'FC'});

%% SCATTERPLOT FC/DFA

 ibif = 14;

figure_white; hold on; title(sprintf('Noise: %d',isig));

for isig = 5 : 5

  scatter(dfa(:,ibif,isig),squeeze(mean(fc(:,:,ibif,isig),1)))

  xlabel('DFA exponent'); ylabel('Mean FC');

  [r,p]=corr(dfa(:,ibif,isig),squeeze(mean(fc(:,:,ibif,isig),1))')

end
%% PLOT STUFF ON CORTEX

dat = zeros(3000,1);
dat2 = zeros(3000,1);
ibi = 6
for i = 1 : 90
  
  load aalmask_grid_cortex3000

  f=fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
  a=textscan(f,'%d %s %d','headerlines',0);

  idx=strcmp(a{2}(i),aalgrid.labels);
  
  dat(idx) = dfa(i,ibi);
  dat2(idx) = dfa(i,ibi+1);

end
  
%%
load sa_meg_template;
dat1 = dat2-dat;


g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;
dd = 0.75;

cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);
%   cmap =   [0.941176474094391 0.941176474094391 0.941176474094391;0.935294091701508 0.935294091701508 0.935294091701508;0.929411768913269 0.929411768913269 0.929411768913269;0.923529386520386 0.923529386520386 0.923529386520386;0.917647063732147 0.917647063732147 0.917647063732147;0.911764681339264 0.911764681339264 0.911764681339264;0.905882358551025 0.905882358551025 0.905882358551025;0.899999976158142 0.899999976158142 0.899999976158142;0.867568910121918 0.852902054786682 0.877736866474152;0.835137844085693 0.805804193019867 0.855473756790161;0.802706778049469 0.758706271648407 0.833210647106171;0.770275771617889 0.711608350276947 0.81094753742218;0.737844705581665 0.664510428905487 0.78868442773819;0.705413639545441 0.617412567138672 0.766421318054199;0.672982573509216 0.570314645767212 0.744158208370209;0.647758424282074 0.533682942390442 0.72684246301651;0.622534275054932 0.497051239013672 0.709526658058167;0.597310066223145 0.460419565439224 0.692210912704468;0.572085916996002 0.423787862062454 0.674895167350769;0.54432600736618 0.383473604917526 0.655838668346405;0.516566097736359 0.343159347772598 0.636782169342041;0.488806158304214 0.30284509062767 0.617725670337677;0.461046248674393 0.262530833482742 0.598669171333313;0.433286309242249 0.222216591238976 0.579612672328949;0.405526399612427 0.181902334094048 0.560556173324585;0.38470646739006 0.151666641235352 0.546263813972473;0.363886535167694 0.121430955827236 0.531971454620361;0.343066573143005 0.0911952629685402 0.51767909526825;0.322246640920639 0.0609595701098442 0.503386676311493;0.301426708698273 0.0307238809764385 0.489094316959381;0.280606776475906 0.000488189951283857 0.474801957607269;0.293446481227875 0.00417250022292137 0.461163550615311;0.306286215782166 0.00785681046545506 0.447525143623352;0.319125920534134 0.0115411207079887 0.433886736631393;0.331965625286102 0.0152254309505224 0.420248329639435;0.344805359840393 0.0189097411930561 0.406609892845154;0.357645064592361 0.0225940514355898 0.392971485853195;0.370484799146652 0.0262783616781235 0.379333078861237;0.383324503898621 0.0299626719206572 0.365694671869278;0.396164208650589 0.0336469821631908 0.352056264877319;0.40900394320488 0.0373312942683697 0.338417857885361;0.421843647956848 0.0410156026482582 0.324779450893402;0.434683352708817 0.044699914753437 0.311141043901443;0.447523087263107 0.0483842231333256 0.297502636909485;0.460362792015076 0.0520685352385044 0.283864229917526;0.473202496767044 0.0557528436183929 0.270225793123245;0.486042231321335 0.0594371557235718 0.256587386131287;0.498881936073303 0.0631214678287506 0.242948994040489;0.511721670627594 0.0668057799339294 0.229310572147369;0.52456134557724 0.0704900845885277 0.215672165155411;0.537401080131531 0.0741743966937065 0.202033758163452;0.550240814685822 0.0778587087988853 0.188395351171494;0.563080549240112 0.0815430209040642 0.174756944179535;0.575920224189758 0.0852273255586624 0.161118522286415;0.588759958744049 0.0889116376638412 0.147480115294456;0.60159969329834 0.0925959497690201 0.133841708302498;0.614439368247986 0.0962802618741989 0.120203301310539;0.627279102802277 0.0999645665287971 0.106564886868;0.640118837356567 0.103648878633976 0.0929264798760414;0.652958512306213 0.107333190739155 0.0792880654335022;0.665798246860504 0.111017502844334 0.0656496584415436;0.678637981414795 0.114701807498932 0.0520112477242947;0.691477656364441 0.118386119604111 0.0383728370070457;0.704317390918732 0.12207043170929 0.024734428152442;0.70988517999649 0.127976059913635 0;0.714490175247192 0.141817703843117 0;0.719095170497894 0.155659362673759 0;0.723700165748596 0.169501006603241 0;0.728305160999298 0.183342665433884 0;0.73291015625 0.197184309363365 0;0.737515151500702 0.211025953292847 0;0.742120146751404 0.224867612123489 0;0.746725142002106 0.238709256052971 0;0.751330137252808 0.252550899982452 0;0.75593513250351 0.266392558813095 0;0.760540127754211 0.280234217643738 0;0.765145123004913 0.294075846672058 0;0.769750118255615 0.307917505502701 0;0.774355113506317 0.321759164333344 0;0.778960108757019 0.335600793361664 0;0.783565163612366 0.349442452192307 0;0.788170158863068 0.363284111022949 0;0.79277515411377 0.377125769853592 0;0.797380149364471 0.390967398881912 0;0.801985144615173 0.404809057712555 0;0.806590139865875 0.418650716543198 0;0.811195135116577 0.432492345571518 0;0.815800130367279 0.446334004402161 0;0.820405125617981 0.460175663232803 0;0.825010120868683 0.474017292261124 0;0.829615116119385 0.487858951091766 0;0.834220111370087 0.501700580120087 0;0.838825106620789 0.515542268753052 0;0.84343010187149 0.529383897781372 0;0.848035097122192 0.543225526809692 0;0.852640092372894 0.557067215442657 0;0.857245087623596 0.570908844470978 0;0.861850082874298 0.584750533103943 0;0.866455078125 0.598592162132263 0;0.871060073375702 0.612433791160583 0;0.875665068626404 0.626275479793549 0;0.880270063877106 0.640117108821869 0;0.884875059127808 0.653958737850189 0;0.88948005437851 0.667800426483154 0;0.894085049629211 0.681642055511475 0;0.898690044879913 0.695483684539795 0;0.903295040130615 0.70932537317276 0;0.907900035381317 0.72316700220108 0;0.912505030632019 0.737008631229401 0;0.917110025882721 0.750850319862366 0;0.921715021133423 0.764691948890686 0;0.926320016384125 0.778533577919006 0;0.930925071239471 0.792375266551971 0;0.935530066490173 0.806216895580292 0;0.940135061740875 0.820058524608612 0;0.944740056991577 0.833900213241577 0;0.949345052242279 0.847741842269897 0;0.953950047492981 0.861583530902863 0;0.958555042743683 0.875425159931183 0;0.963160037994385 0.889266788959503 0;0.967765033245087 0.903108477592468 0;0.972370028495789 0.916950106620789 0;0.97697502374649 0.930791735649109 0;0.981580018997192 0.944633424282074 0;0.986185014247894 0.958475053310394 0;0.990790009498596 0.972316682338715 0;0.995395004749298 0.98615837097168 0;1 1 0];
par_interp = spatfiltergauss(dat1,g1,dd,g2);

figure; set(gcf,'color','white'); hold on;

  para = [] ;

  para.colorlimits = [-.02 .02];
  
  for iplot = 1 : 6
      
    subplot(3,2,iplot)

     para.myviewdir = viewdir(iplot,:);
     a = sa_meg_template.cortex10K;
     
     if iplot == 5
      	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
     elseif iplot == 6
      	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
     end
     
    pconn_showsurface(a,para,par_interp)
    colormap(cmap)
    camlight headlight

  end




%%

struct(:,:,:,1) = coh;
struct(:,:,:,2) = fc_sim;
struct(:,:,:,3) = repmat(SC,[1 1 21]);

save('~/ampcorr_coh_sc.mat','struct');

%%


figure; set(gcf,'color','white'); hold on

for i = 1 : 45
  
  a = fc_sim(:,:,i);
  b = coh(:,:,i);
  
  r(i) = corr(a(:),b(:));
  
  ra(i) = corr(a(:),SC(:));
  rb(i) = corr(b(:),SC(:));
  
  m(i) = nanmean(a(:));
  n(i) = nanmean(b(:));
  
end

figure; set(gcf,'color','w'); hold on


plot(para.sig(1:45),ra,'r','linewidth',3)
plot(para.sig(1:45),rb,'b','linewidth',3)

xlabel('Noise');
ylabel('Correlation(amp. corr,abs(coh))');


figure; set(gcf,'color','w'); hold on

plot(para.sig(1:45),r,'b','linewidth',3)

xlabel('Noise');
ylabel('Correlation(amp. corr,abs(coh))');

figure; set(gcf,'color','w'); hold on
plot(m,'r','linewidth',3)
plot(n,'b','linewidth',3)

xlabel('Noise');
ylabel('Mean connectivity')


%%
figure; set(gcf,'color','white');

% PLOT CONNECTIVITY MATRICES

subplot(2,3,1);

imagesc(coh);  axis square;

subplot(2,3,2);

imagesc(fc_sim);  axis square;

subplot(2,3,3);

icoh = imag([(cs).^2]./[diag(cs)*diag(cs)']);

imagesc(icoh);  axis square;

% PLOT CORRELATIONS

subplot(2,3,4);

scatter(coh(:),fc_sim(:)); axis square;  box on;
r=corr(coh(:),fc_sim(:));

subplot(2,3,5);

scatter(coh(:),icoh(:)); axis square;  box on;
r=corr(coh(:),icoh(:));

subplot(2,3,6);

scatter(fc_sim(:),icoh(:)); axis square; box on;
r=corr(fc_sim(:),icoh(:));


fprintf('Amp. corr with SC: %.2f\n',corr(fc_sim(:),SC(:)));

fprintf('Coh. with SC: %.2f\n',corr(coh(:),SC(:)));

% function coh = comp_coherence(f,nchan,nseg)
%
% nchan = 90;
% f = 11;
% seglen  = 400; % in samples
% nseg    = floor(size(xs,1)/seglen);
% overlap = 0.5;
% win = hanning(seglen);
% cs=zeros(nchan,nchan,nseg);
%
% for iseg = 1 : nseg
%
%   idx = (iseg-1)*seglen*overlap+1:(iseg-1)*seglen*overlap+seglen;
%
%   dat = xs(idx,:);
%
%   datfft = fft(repmat(win,[1 size(dat,2)]).*dat);
%
%   cs(:,:,iseg)      = cs(:,:,iseg)+conj(datfft(f,:)'*datfft(f,:));
%   %   abscoh(:,:,iseg)  = cs(:,:,iseg)./sqrt(diag(cs(:,:,iseg))*diag(cs(:,:,iseg))');
%
%
% end
%
% cs = nanmean(cs,3);
%
% coh(:,:,isig) = [abs(cs).^2]./[abs(diag(cs))*abs(diag(cs))'];
%


