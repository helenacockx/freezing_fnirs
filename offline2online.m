function [data_combi, events]=offline2online(filename_online, varargin)

%% initialisation
vis=ft_getopt(varargin, 'vis', true);
delay_devices=ft_getopt(varargin, 'delay_devices'); %24065-24068 in samples

[path, name, ext]=fileparts(filename_online);
fparts=regexp(name, 'sub-(?<sub>\w+)_task-gait_acq-(?<acq>\w+)_rec-(?<rec>\d+)_nirs', 'names');
filename_24065 = fullfile(path, sprintf('sub-%s_task-gait_acq-24065_rec-%s_nirs.oxy3', fparts.sub, fparts.rec));
filename_24068 = fullfile(path, sprintf('sub-%s_task-gait_acq-24068_rec-%s_nirs.oxy3', fparts.sub, fparts.rec));

% %%  initialisation 
% subID='HC57';
% subject= sprintf('sub-%s', subID);
% rec = 1;
% root_dir='C:\Users\helen\Documents\freezing_fnirs\data\';
% 
% filename_24065 = fullfile(root_dir, 'source_standard', subject, 'nirs', sprintf('%s_task-gait_acq-24065_rec-%02d_nirs.oxy3', subject, rec));
% filename_24068 = fullfile(root_dir, 'source_standard', subject, 'nirs', sprintf('%s_task-gait_acq-24068_rec-%02d_nirs.oxy3', subject, rec));
% filename_online = fullfile(root_dir, 'source_standard', subject, 'nirs',sprintf('%s_task-gait_acq-online_rec-%02d_nirs.oxy4', subject, rec));
% filename_output = fullfile(root_dir, 'raw', subject, 'nirs', sprintf('%s_task-gait_rec-%02d_nirs.snirf', subject, rec));
% 
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\oxysoft2matlab-v1.69\priv
% addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
% cd(fullfile(root_dir, 'source_standard', subject, 'nirs'))

%% read data
% online data
cfg=[];
cfg.dataset=filename_online;
data_online=ft_preprocessing(cfg);

% offline data 24065 (==> gather metainfo form readoxy3?)
[rawOD_24065,metaInfo_24065,ADvalues_24065] = readoxy3file(filename_24065);
global ADC; % ADvalues of readoxy3file.m are converted in line 113
ADvalues_24065=ADC'; 
metaInfo_24065.OptodeTemplateID = 131;
metaInfo_24065.DPF = [6 6]; % see online file
metaInfo_24065.Position = [30 30]; % see optodetemplates.xml (DistanceScale)
metaInfo_24065.abs_K = [1.1 1.1]; % not sure what this means, taken over from online file
metaInfo_24065.abs_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24065.use_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24065.Gradient = [5 5];% not sure what this means, taken over from online file
[oxyOD2_24065, metaInfo_24065] = arrangeoxy3optodes(rawOD_24065', metaInfo_24065);

% offline data 24068
[rawOD_24068,metaInfo_24068,ADvalues_24068] = readoxy3file(filename_24068);
global ADC;
ADvalues_24068=ADC';
metaInfo_24068.OptodeTemplateID = 130;
metaInfo_24068.DPF = [6 6]; % see online file
metaInfo_24068.Position = [30 30]; % see optodetemplates.xml (DistanceScale)
metaInfo_24068.abs_K = [1.1 1.1]; % not sure what this means, taken over from online file
metaInfo_24068.abs_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24068.use_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24068.Gradient = [5 5];% not sure what this means, taken over from online file
[oxyOD2_24068, metaInfo_24068] = arrangeoxy3optodes(rawOD_24068', metaInfo_24068);

%% map offline channels onto online channels
nsamples=min(metaInfo_24065.nbSamples, metaInfo_24068.nbSamples);
nchan=sum(strcmp(data_online.hdr.chantype, 'nirs'));
OD_offline=nan(nchan, nsamples);
% 24065
for k=1:numel(metaInfo_24065.Sys.subtemplate) % loop over subtemplates
  for i=1:length(metaInfo_24065.Sys.subtemplate(k).RxTx)
    idx=find(contains(data_online.label, metaInfo_24065.Sys.subtemplate(k).RxTx{i}));
    OD_offline(idx,:)=oxyOD2_24065.Sys.subtemplate(k).OD(1:nsamples,2*i-1:2*i)';
  end
end
% 24068
for k=1:numel(metaInfo_24068.Sys.subtemplate) % loop over subtemplates
  for i=1:length(metaInfo_24068.Sys.subtemplate(k).RxTx) % now Rx starts from Rx9 and Tx from Tx11...
    try
      parts=regexp(metaInfo_24068.Sys.subtemplate(k).RxTx{i}, 'Rx(?<RxID>\d+)-Tx(?<TxID>\d+)(?<Txsub>\D)', 'names');
      chan=sprintf('Rx%d-Tx%d%s', str2double(parts.RxID)+8, str2double(parts.TxID)+10, parts.Txsub);
    catch
      parts=regexp(metaInfo_24068.Sys.subtemplate(k).RxTx{i}, 'Rx(?<RxID>\d+)-Tx(?<TxID>\d+)', 'names');
      chan=sprintf('Rx%d-Tx%d', str2double(parts.RxID)+8, str2double(parts.TxID)+10);
    end
    idx=find(contains(data_online.label, chan));
    OD_offline(idx,:)=oxyOD2_24068.Sys.subtemplate(k).OD(1:nsamples,2*i-1:2*i)';
  end
end
% remark: this is in fact the same order as in the data_online file...

%% find offset of offline vs online files and 24065 vs 24068
% online files
event_online{1}=trig_detect(data_online.trial{1}, 1, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 1 (24065)
event_online{2}=trig_detect(data_online.trial{1}, 53, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 53 (24068)
for k=1:2 % find for both devices
    for i=1:length(event_online{k})-1
      if contains(event_online{k}(i).type, 'up') && contains(event_online{k}(i+1).type, 'down') && (event_online{k}(i+1).sample-event_online{k}(i).sample>50 & event_online{k}(i+1).sample-event_online{k}(i).sample<1000)
        fprintf('upgoing flank detected at sample %d, verify in the figure if this is correct \n', event_online{k}(i).sample)
        break
      end
    end
  trig_online{k}=event_online{k}(i).sample;
end
if vis
  figure; subplot(2,1,1); plot(data_online.trial{1}([1 53],1:15000)'); hold on; plot(trig_online{1}, 4, 'o'); plot(trig_online{2}, 4, 'o'); title('trigger detection online data'); legend({'24065', '24068'});
end

event_offline{1}=trig_detect(OD_offline, 1, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 1 (24065)
event_offline{2}=trig_detect(OD_offline, 53, 1, 15000, 3);% find up and down going flank in the first 5 minutes in channel 53 (24068)
for k=1:2 % find for both devices
  for i=1:length(event_offline{k})-1
    if contains(event_offline{k}(i).type, 'up') && contains(event_offline{k}(i+1).type, 'down') && (event_offline{k}(i+1).sample-event_offline{k}(i).sample>50 & event_offline{k}(i+1).sample-event_offline{k}(i).sample<1000)
      fprintf('upgoing flank detected at sample %d, verify in the figure if this is correct \n', event_offline{k}(i).sample)
      break
    end
  end
  trig_offline{k}=event_offline{k}(i).sample;
end
if vis
  subplot(2,1,2); plot(OD_offline([1 53],1:15000)'); hold on; plot(trig_offline{1}, 4, 'o'); plot(trig_offline{2}, 4, 'o'); title('trigger detection offline data'); legend({'24065', '24068'});
end

% calculate offsets
offset_24065=trig_offline{1}-trig_online{1};
offset_24068=trig_offline{2}-trig_online{2};

% correct for shift between devices
if delay_devices % = 24068-24065
  fprintf('correcting for shift between devices with offset of %d samples', delay_devices)
  offset_24065=offset_24065-delay_devices;
end

% % offsets between 24065 vs 24068
% if ~isempty(metaInfo_24065.Event) & ~isempty(metaInfo_24068.Event)
%   try
%     offset_devices=metaInfo_24065.Event{1}-metaInfo_24068.Event{1}; % if multiple events
%   catch
%     offset_devices=metaInfo_24065.Event-metaInfo_24068.Event;
%   end
%   fprintf('Correcting offset of 24065 with %d samples based on the synchronisation event that was inserted in the two offline devices (offset between 24065-24068 = %d samples) \n', offset_24065-offset_24068-offset_devices, offset_devices)
%   offset_24065=offset_24068+offset_devices;
% else
%   fprintf('No synchronisation event was inserted in the offline files. Please check if both devices are well synchronized \n');
% end  
% offsets between 24065 vs 24068

if isfield(metaInfo_24065, 'Event') & isfield(metaInfo_24068, 'Event') & ~strcmp(fparts.sub, 'HC42') % erroneuos events were inserted in HC42
  try
    if length([metaInfo_24065.Event{:}])==length([metaInfo_24068.Event{:}])
      delays_offline=[metaInfo_24068.Event{:}]'-[metaInfo_24065.Event{:}]';
      table(delays_offline)
    end
    if strcmp(fparts.sub, 'HC76') % errnuous first 24065 event
       offset_devices=metaInfo_24068.Event{1}-metaInfo_24065.Event{2};
    else 
      offset_devices=metaInfo_24068.Event{1}-metaInfo_24065.Event{1}; % if multiple events
    end
  catch
    offset_devices=metaInfo_24068.Event-metaInfo_24065.Event;
  end
  if abs((offset_24068-offset_24065)-offset_devices)>25 % if devices are off with more than 0.5 seconds, throw error
    warning('Desynchronization between the devices with more than 0.5 seconds')
  else
    fprintf('Based on the offline events, there might be a delay between the devices in the online file of %.03f seconds (delay in offline file was %d samples). \n', (abs(offset_devices-(offset_24068-offset_24065)))/50, offset_devices);
  end
  fprintf('correcting for this delay... \n')
  delay_devices=offset_devices-(offset_24068-offset_24065);
  offset_24068=offset_24068+delay_devices;
else
  figure; plot(ADvalues_24065(3, 50:end)); hold on; plot(ADvalues_24068(3, 50:end)); title(sprintf('sub-%s rec-%s: 24065 vs 24068', fparts.sub, fparts.rec)); legend({'24065', '24068'});
  fprintf('No synchronisation event was inserted in the offline files. Please check in the figure what the delay is between the two offline files \n');
  pause;  confirm=0;
  while confirm~=1
    offset_devices=input('What is the delay between 24068-24065 (in samples)? \n')
    confirm = input('Is this the definite answer? 1/0 \n');
  end
  delay_devices=offset_devices-(offset_24068-offset_24065);
  offset_24068=offset_24068-delay_devices;
end  

% apply offsets to the offline data 
fprintf('Applying offset of %d samples for offline file 24065 and of %d samples for offline file 24068 \n', offset_24065, offset_24068)
nsamples=size(OD_offline,2)-max(offset_24065, offset_24068);
data=nan(length(data_online.label),nsamples);
data(1:52,:)=OD_offline(1:52, offset_24065+1:offset_24065+nsamples); % 24065
data(53:96,:)=OD_offline(53:96, offset_24068+1:offset_24068+nsamples); % 24068
data(97:100,:)=ADvalues_24065(:,offset_24065+1:offset_24065+nsamples);
data(101:104,:)=ADvalues_24068(:,offset_24068+1:offset_24068+nsamples);
% and paste into online_data
data_combi=data_online;
data_combi.sampleinfo=[1 nsamples];
data_combi.time{1}=[0:1/data_combi.fsample:(nsamples-1)/data_combi.fsample];
data_combi.trial{1}=data;
data_combi.hdr.nSamples=nsamples;

% check data
if vis
  figure; subplot(3,1,1); plot(data_online.trial{1}(3,50:end), 'b-'); hold on; plot(data_combi.trial{1}(3,50:end), 'g.'); title(sprintf('sub-%s rec-%s: 24065 vs online', fparts.sub, fparts.rec)); legend({'offline data', 'online data'});
  subplot(3,1,2); plot(data_online.trial{1}(53,50:end), 'b-'); hold on; plot(data_combi.trial{1}(53,50:end), 'g.'); title(sprintf('sub-%s rec-%s: 24068 vs online', fparts.sub, fparts.rec)); legend({'offline data', 'online data'});
  subplot(3,1,3); plot(data_combi.trial{1}(99, 50:end)); hold on; plot(data_combi.trial{1}(103, 50:end)); title(sprintf('sub-%s rec-%s: 24065 vs 24068', fparts.sub, fparts.rec)); legend({'24065', '24068'});
end

%% read in events
events=ft_read_event(filename_online, 'chanindx', -1, 'type', 'event');
% check with offline events if existent
start_runs=events(strcmp({events.value}, 'LSL start_run'));
try
if length(start_runs)==length(metaInfo_24068.Event)
  delays_event=([metaInfo_24065.Event{:}]'-offset_24065)-[start_runs.sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event{:}]'-offset_24068)-[start_runs.sample]';
  table(delays_event)
  
  delays_event=([metaInfo_24065.Event{1}]'-offset_24065)-[start_runs(1).sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event{1}]'-offset_24068)-[start_runs(1).sample]';
  table(delays_event)
  
    delays_event=([metaInfo_24065.Event]'-offset_24065)-[start_runs(1).sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event]'-offset_24068)-[start_runs(1).sample]';
  table(delays_event)
  
    delays_event=([metaInfo_24065.Event{[1:3 5]}]'-offset_24065)-[start_runs([1:3 6]).sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event{:}]'-offset_24068)-[start_runs([1:3 6]).sample]';
  table(delays_event)
  
      delays_event=([metaInfo_24065.Event{[2:4]}]'-offset_24065)-[start_runs([2:4]).sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event{:}]'-offset_24068)-[start_runs([2:4]).sample]';
  table(delays_event)
  
        delays_event=([metaInfo_24065.Event{:}]'-offset_24065)-[start_runs([1 3:4]).sample]';
  table(delays_event)
  delays_event=([metaInfo_24068.Event{:}]'-offset_24068)-[start_runs([1 3:4]).sample]';
  table(delays_event)
end
end

%% check data
if vis
  cfg                = [];
  cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
  cfg.viewmode       = 'vertical';
  cfg.channel = 'Rx*';
  cfg.event =events;
  cfg.ploteventlabels= 'colorvalue';
  cfg.plotlabels= 'yes';
  cfg.fontsize=5;
  cfg.continuous     = 'yes';
  cfg.blocksize  = 60;
  cfg.nirsscale =100;
  ft_databrowser(cfg, data_combi);
end


%% SUBFUNCTIONS %%
function event=trig_detect(dat, chan, beginsmp, endsmp, threshold)

% select data
dat=dat(chan, beginsmp:endsmp);
% discretize the signal
dat(dat< threshold) = 0;
dat(dat>=threshold) = 1;
% convert the trigger into an event with a value at a specific sample
event=[];
begpad = dat(1);
endpad = dat(end);
difftrace = diff([begpad dat endpad]);

for j=find(difftrace~=0)
  if difftrace(j)>0
    event(end+1).type   = 'up';        % distinguish between up and down flank
    event(end  ).sample = j + beginsmp - 1;      % assign the sample at which the trigger has gone up
    event(end  ).value  = dat(j);      % assign the trigger value just _after_ going up
  elseif difftrace(j)<0
    event(end+1).type   = 'down';      % distinguish between up and down flank
    event(end  ).sample = j + beginsmp - 1;      % assign the sample at which the trigger has gone down
    event(end  ).value  = dat(j-1);    % assign the trigger value just _before_ going down
  end
end
