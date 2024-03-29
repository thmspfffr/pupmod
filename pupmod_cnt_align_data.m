%% pupmod_cnt_align_data
% subj < 25: all good
% difficulty with subj > 25: i first rejected artifacts, and then removed
% first 20s of the data. 

clear

v = 1;


addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

% SUBJECT10, M3, B1: no triggers
% sth was weird with subject 24, also subject 10 failed

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = SUBJLIST
  for m = 1:3
    for iblock = 1:2

     	load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v1.mat',isubj,m,iblock))
      if isubj > 24
%         data.trial{1} = [nan(size(data.trial{1},1),8000) data.trial{1}];
      end
      
      try
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_trig_raw_s%d_m%d_b%d_v1.mat',isubj,m,iblock))
      catch me
        warning(sprintf('No triggers for s%d, m%d, b%d',isubj,m,iblock))
      end
                
      if isubj <= 24 
        load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_artvec_s%d_m%d_b%d_v1.mat',isubj,m,iblock));
        % do this, because artvec for s8 is somehow incorrect. data.idx
        % provides correct artifacts. if art vec is reduced like below,
        % everythign is okay
        if isubj == 8 && m == 2 && iblock == 1
          art(5,:)=[];
        elseif isubj == 8 && m == 2 && iblock == 2
          art(4:6,:)=[];      
        end
      else
        load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_artvec_s%d_m%d_b%d_v1.mat',isubj,m,iblock));
      end
      
      %	--------------------------
      % START REALIGNING MEG SIGNAL
      %	--------------------------
      art = round(art./3);
      art(art==0)=1;
      
      if isubj==10 && m == 3
        trig.trial{1} = zeros(1,trig.sampleinfo(2));
        trig.trial{1} = resample(trig.trial{1},400,1200);
      else
        cfg = [];
        cfg.resamplefs = 400;
        trig = ft_resampledata(cfg,trig);
      end

      artvec = zeros(size(trig.trial{1}));
      for iart = 1 : size(art,1)
        artvec(art(iart,1):art(iart,2))=1;
      end
      bw = bwlabel(artvec);
      difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
      
      fprintf('Difference in data: %d samples\n',difference)
      
      if (size(data.trial{1},2)+sum(bw>0))~=size(data.idx(:),1)
        difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
        fprintf('Diff = %d | Art = %d\n',difference,length(art))
        if difference > 0 && difference <= 25      
          while difference > 0 
            i = randi(max(bw));
            bw(find(bw==i,1,'first'))=0;
            difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
          end
        elseif difference < 0
          while difference > 0 
            i = randi(max(bw));
            bw(find(bw==i,1,'last')+1)=i;
            difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
          end
        else
          warning('Many artifacts? Differs by more than 25 samples')
          while difference > 0 
            i = randi(max(bw));
            bw(find(bw==i,1,'first'))=0;
            difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
          end
        end
      end
   
      % --------------------------
      % FIND TRIGGER IN DATA
      % --------------------------
      trig_idx = bwlabel(abs(trig.trial{1})>0.5);
      
      if round((round(mean(find(trig_idx==max(trig_idx))))-round(mean(find(trig_idx==1))))/trig.fsample)~=600
        start = round(mean(find(trig_idx==1)));
        off = start + 600*trig.fsample;
        fprintf('S%dM%dB%d: Only start trigger\n',isubj,m,iblock)
      elseif round((round(mean(find(trig_idx==max(trig_idx))))-round(mean(find(trig_idx==1))))/trig.fsample)==600
        fprintf('S%dM%dB%d: Both triggers\n',isubj,m,iblock)
        start = round(mean(find(trig_idx==1)));
        off = round(mean(find(trig_idx==max(trig_idx))));
      else
        warning('error')
      end
      
      % --------------------------
      % CREATE NEW DATA SET, WITH NANS WHERE ARTIFACTS WERE
      % --------------------------
      clear  dat_len art_len
       
      for ibw = 1 : max(bw)
        if ibw == 1
          dat_len(1,:) = length(1:find(bw==ibw,1,'first')-1);
        else
          dat_len(ibw,:) = length(find(bw==ibw-1,1,'last')+1:find(bw==ibw,1,'first')-1);
        end
        art_len(ibw,:) = sum(bw==ibw);
      end
      dat_len(ibw+1,:) = length(find(bw==ibw,1,'last')+1:length(bw));
      if dat_len(ibw+1,:) == 0 
        dat_len(ibw+1,:) = [];
      end
      
      dat_new = [];
      data_old = data.trial{1}; data.trial = [];
      
      for i = 1 : size(dat_len,1)
        if i <= size(art_len,1)
          dat1 = data_old(:,1:dat_len(i));
          data_old(:,1:dat_len(i)) = [];
          dat_new = [dat_new dat1 nan(size(dat1,1),art_len(i))];
        else
          dat_new = [dat_new data_old];
        end
      end
      
      data.trigger = trig.trial{1};
      data.art = artvec;
      
      data.start_of_recording = start;
      data.end_of_recording = off;
      
      data.time =  1/data.fsample:1/data.fsample:size(dat_new,2)/data.fsample;
      data.trial = dat_new;
      data.cfg = [];
      
      save(sprintf('~/pupmod/proc/pupmod_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v),'data')
      
      clear data dat_new data_old dat1 start off bw art trig data artvec  
      
    end
  end
end
    
  
%% LOAD DATA 

for isubj = SUBJLIST
  for m = 1:3
    for iblock = 1:2
      
       fprintf('Processing s%d m%d b%d ...\n', isubj,m,iblock)

      try
        load(sprintf('~/pupmod/proc/pupmod_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v))
      catch me
        continue
      end
      
      if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
        if (data.end_of_recording-600*data.fsample)<1
          data.start_of_recording = 1;
        else
          data.start_of_recording = data.end_of_recording-600*data.fsample;
        end
        warning('Only End trigger')
      elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
        if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
          data.end_of_recording = size(data.trial,2);
        else
          data.end_of_recording = data.start_of_recording+600*data.fsample;
        end
        warning('Only Start trigger')
      elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
        warning('No triggers')
      elseif isnan(data.start_of_recording) && isnan(data.end_of_recording)
        data.start_of_recording = 5000;
        data.end_of_recording = 235000;
        warning('NaN triggers')
      end
      
      dat = data.trial(:,data.start_of_recording:data.end_of_recording);
      
      save(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v),'dat');

    end
  end
end
    
