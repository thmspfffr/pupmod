%% pupmod_src_hopf
% pupmod_src_hopf_fit.m

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 1;
fsample  = 400;
para.bif = -0.15:0.01:0.15;
para.G   = 0.2:0.2:5;
para.sig = 0.01;
para.freq = 10;
% --------------------------------------------------------

para.freq = para.freq(para.freq<40);

outdir   = '/home/tpfeffer/pupmod/proc/';

%%
clear r n
v = 1;
cnt = 0;

for isig = 1  : 1%length(para.sig)
  for ifreq = 1 : length(para.freq);
    for ig = 1 : length(para.G)
      ig
      for ibif = 1 : length(para.bif)
%         try
          load(sprintf('/home/tpfeffer/pupmod/proc/pupmod_fc_sim_f%d_s%d_b%d_g%d_v%d.mat',ifreq,isig,ibif,ig,v))
          
          m(isig,ifreq,ig,ibif) = mean(nonzeros(triu(nanmean(fc_sim,3),1)));
          
%           [n(:,ibif,ig,isig),k]=hist(nonzeros(triu(fcd,1)),-1:0.01:1);
          
      end
    end
  end
%   save([outdir sprintf('pupmod_src_hopf_simcomb_s%d_v%d.mat',isig,v)],'fc_all');
  disp('done');
end

fprintf('Done loading. Cound: %d',cnt)
error('!')
%%
isig = 1;
v_corr = 8;
ord = pconn_randomization;

load([outdir sprintf('pupmod_src_hopf_simcomb_s%d_v%d.mat',isig,v)]);

for isubj = [4 5 6 7 9 10 11 12 13 15 16 19 21 23 24]
  for iblock = 1 : 2
    for m = 1 : 3
      im = find(ord(isubj,:)==m);
%     figure_white
      for ifreq = 11 : 1
   
        
        %%
        ifreq =  6
       
      load(sprintf([outdir 'pupmod_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifreq,v_corr));

      fprintf('Computing FCD for s%d b%d m%d f%d ...\n',isubj,iblock,m,ifreq)
      for iep = 1 : 100%size(powcorr,3)
        for jep = 1 : 100%size(powcorr,3) 
            
          fcd_emp(iep,jep) = corr(nonzeros(triu(powcorr(:,:,iep),1)),nonzeros(triu(powcorr(:,:,jep),1)));
          
        end
      end
    	fprintf('Computing FCD for s%d b%d m%d f%d ... Done!\n',isubj,iblock,m,ifreq)
imagesc(fcd_emp,[-0.2 .2]); colormap(cmap)
%%

      clear r d
      for ig = 4% : length(para.G)
        fprintf('Computing freq%d, G%d ...\n',ifreq,ig)
        for ibif = 1 : length(para.bif)
                    
          fc_sim =  triu(fc_all(:,:,ibif,ig),1);
          
          r(ibif) = corr(tmp(:),fc_sim(:));
          d(ibif) = mean(abs(tmp(:)-fc_sim(:)));
          
          clear fc_sim
          
        end
      end
      
      [~, mm(isubj,iblock,ifreq,m)] = max(r(1:16));
      
    end
      
%       subplot(4,4,ifreq)
      
%       imagesc(r,[-0.75 0.75]); colormap(jet);
%       
% %       set(gca,'xtick',[1:5:length(para.bif)],'xticklabel',para.bif(1:5:length(para.bif)));
%             set(gca,'xtick',[]);
% 
% %       set(gca,'ytick',[1:5:length(para.G)],'yticklabel',para.G(1:5:length(para.G)));
%             set(gca,'ytick',[]);
% 
%       set(gca,'tickdir','out')
%       
%       xlabel('Bifurcation parameter'); ylabel('Global coupling');
      
%       print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_hopf_parameters_s%d_b%d_f%d_v%d.jpg',isubj,iblock,para.freq,v))
    end
  end
end