%% pupmod_src_hopf_gain

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 1;
fsample  = 400;
para.bif = -0.2:0.02:0.04;
para.G   = 0:2:14;
para.sig = 0.01;
para.freq = 10;
para.gain     = 0.02:0.02:0.2;
% --------------------------------------------------------

para.epleng = 60;
para.overlap = 0.98;

para.freq = 40*ones(90,1);%para.freq(para.freq<40);

outdir   = '/home/tpfeffer/pupmod/proc/';

%% RUN HOPF MODEL
for ifoi = 1 : 1
  for igain = 1 : length(para.gain)
  for isig = 1 : length(para.sig)
    for ibi = 1 : length(para.bif)
      for iG = 1 : length(para.G)
        
        if ~exist(sprintf([outdir 'pupmod_src_hopf_gain_f%d_s%d_b%d_g%d_v%d_processing.txt'],ifoi,isig,ibi,iG,v))
          system(['touch ' outdir sprintf('pupmod_src_hopf_gain_f%d_s%d_b%d_g%d_v%d_processing.txt',ifoi,isig,ibi,iG,v)]);
        else
          continue
        end
        
        clear xs kura
        
        fs = 400;
        
        load ~/Downloads/SC_AAL_90.mat
        
        areas       = 90;
        C           = SC/max(max(SC))*0.2;
        
        time_steps  = 600;
        
        a           = para.bif(ibi);  % bifurcation param.
        we          = para.G(iG);     % coupling
        sig         = para.sig(isig);           % amp. of noise?
        f           = para.freq;             % freq. = omega
        
        dt          = 1/fsample;
        
        amp         = zeros(time_steps/dt,areas);
        %     ys          = zeros(time_steps/dt,areas);
        
        kura        = zeros(time_steps/dt,1);
        
        x           = 0.1 * randn(areas,1);
        y           = 0.1 * randn(areas,1);
        nn          = 0;
        omega       = (2*pi).*f;
        
        gain        = para.gain(igain);
        
        Isubdiag = find(tril(ones(areas),-1));
        
        fprintf('Initializing...\n');
        
        for t=0:dt:100
          
          sumax = C*x-diag(C*repmat(x',areas,1));
          sumay = C*y-diag(C*repmat(y',areas,1));
          
          % *supercritical*
          x = x+dt*(resp_fun((a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1),gain)-0.5);
          y = y+dt*(resp_fun((a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1),gain)-0.5);
          
        end
        
        fprintf('Running model...\n');
        
        for t=0:dt:time_steps  %32000
          
          fprintf('Time step %.3f / %d ...\n',t,time_steps);
          
          sumax = C*x-diag(C*repmat(x',areas,1));
          sumay = C*y-diag(C*repmat(y',areas,1));
          
          x = x+dt*(resp_fun((a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1),gain)-0.5);
          y = y+dt*(resp_fun((a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1),gain)-0.5);
          
          %       if mod(t,2)==0
          nn=nn+1;
          %       xs(nn,:)=x';
          %       ys(nn,:)=y';
          amp(nn,:) = sqrt(x.*x+y.*y);
          ri=amp(nn,:)';
          rimax=max(ri);
          ri=ri/rimax;
          kura(nn)=abs(sum(ri.*complex(cos(angle(complex(x,y))),sin(angle(complex(x,y)))))/areas);
          %       end
          
        end
        
        epleng = fsample*para.epleng;
        epshift = epleng-(para.overlap*epleng);
        nep =floor((fsample*time_steps-epleng)/epshift+1);
        
        for iep = 1 : nep
        
          fc_sim(:,:,iep) = corrcoef(amp((iep-1)*epshift+1:(iep-1)*epshift+epleng,:));
          
        end
        
        for iep = 1 : nep
          for jep = 1 : nep

            fcd(iep,jep) = corr(nonzeros(triu(fc_sim(:,:,iep),1)),nonzeros(triu(fc_sim(:,:,jep),1)));

          end
        end
        
        save(sprintf([outdir '/pupmod_fc_sim_gain_f%d_s%d_b%d_g%d_gain%d_v%d.mat'],ifoi,isig,ibi,iG,igain,v),'fcd','fc_sim','kura');
        
      end
    end
  end
end
end
error('!')

%% FCD

v = 2;

for ifoi = 1 : 1
  for igain = 1 : length(para.gain)
    for isig = 1 : length(para.sig)
      for ibi = 1 : length(para.bif)
        for iG = 1 : length(para.G)
          
          load(sprintf([outdir '/pupmod_fc_sim_gain_f%d_s%d_b%d_g%d_gain%d_v%d.mat'],ifoi,isig,ibi,iG,igain,v))
          
          m(isig,ifreq,ig,ibif,igain) = mean(nonzeros(triu(nanmean(fc_sim,3),1)));

        end
      end
    end
  end
end

%%
v_corr = 10;
ord = pconn_randomization;

pc = zeros(90,90,24,3,2);
SUBJLIST =  [4 5 6 7 9 10 11 12 13 15 16 19 21 23 24];
b = inf; clear r;

for isubj = 4%SUBJLIST
  isubj
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      for ifreq = 6 : 6
        
        load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifreq,v_corr));
        emp = nonzeros(triu(nanmean(powcorr,3)));
        
        for ibif =1 : 30
%           ibif
          for ig = 1 : 15%length(para.G)
            isig = 2; ifoi = 9;
            load(sprintf([outdir '/pupmod_fc_sim_f%d_s%d_b%d_g%d_v%d.mat'],ifoi,isig,ibif,ig,v));
            
            fc_sim(fc_sim==1)=0;
            sim = nonzeros(triu(nanmean(fc_sim,3)));
            

            r(ibif,ig,isubj,m,iblock,1)=sum(abs(sim-emp));
            r(ibif,ig,isubj,m,iblock,2)=corr(sim,emp);

            if abs(corr(sim,emp))<abs(b)
              b=corr(sim,emp)   
            end
             
              
          end
          
          
        end
      end
    end
  end
end

nanmean(r(:,:,SUBJLIST,:,:),5)

pc = nanmean(pc(:,:,SUBJLIST,:,:),5);

%%
% BIF = 2;
% 
% clear m
% for isig = 3 
%   
%       for ifreq = 9 : 9
%         
%         
%         r(ibif,ig)=corr(t(:),tmp(:));
% %         figure;
%         m(ibif,ig) = mean(fc_sim(:));
% %         imagesc(fc_sim,[0 1]); drawnow
%         
%       end
%     end
%   end
% %   [~,idx(isig)]=max(m); clear m
% end


        
