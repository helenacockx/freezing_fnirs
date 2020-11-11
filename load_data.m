%% initialisation
close all; clear all
ID='PD11';
sub=['sub-' ID];
root_dir='C:\Users\helen\Documents\freezing_fnirs\data';
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
addpath(fullfile(root_dir, 'scripts'))
dbstop if error
% completed subjects: PD11, HC66

runs=[1:5];
num_run=5;
num_block=13; 

%nirs
delay_devices=0;

% make processed folder
subfolders={'nirs', 'motion', 'video', 'stim'};
for i=1:length(subfolders)
  mkdir(fullfile(root_dir, 'processed', sub, subfolders{i}));
end
% log file
logname=fullfile(root_dir, 'processed', sub, sprintf('%s_load_data.log', sub)); diary(logname);

%% load REC variable
fprintf('\n.........loading lsl data.........\n')
triggerinfo=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*rec-0*_triggerinfo.mat'));
load(fullfile(triggerinfo.folder, triggerinfo.name)); 

% check number of events
lsl_streams={'lsldert0_events', 'lsldert1_events', 'lsldert2_events', 'lsldert3_events', 'lsldert4_events'};
sync=nan(1,length(lsl_streams));
for i=1:length(lsl_streams)
  if ~isequal(sum(contains(REC.(lsl_streams{i}).Data,'TTL_on')),sum(contains(REC.(lsl_streams{i}).Data,'TTL_off')),sum(contains(REC.(lsl_streams{i}).Data,'start_run')),sum(contains(REC.(lsl_streams{i}).Data,'stop_run')), num_run)
    warning('The number of start_runs/TTL_on and stop_runs/TTL_off of %s were not as expected', lsl_streams{i})
  end
  if ~isequal(sum(contains(REC.(lsl_streams{i}).Data,'start_block')),sum(contains(REC.(lsl_streams{i}).Data,'stop_block')), num_block)
    warning('The number of start_blocks and stop_blocks of %s were not as expected', lsl_streaminfos{i})
  end
  sync(i)=sum(contains(REC.(lsl_streams{i}).Data, 'sync'));
end
if 
if length(unique(sync))~=1
  warning('The number of sync events were not equal between the lsl streams')
end
% check delay between events and calculate corrections
lsl_events=test_delay(REC);


%% load nirs data
fprintf('\n.........loading nirs data.........\n')
oxy4files=dir(fullfile(root_dir, 'source_standard', sub, 'nirs', '*acq-online*.oxy4'));
[data_combi, nirsevents]=offline2online(fullfile(oxy4files.folder, oxy4files.name), 'delay_devices', delay_devices);
% temporary save data_combi
fparts=regexp(oxy4files.name, 'sub-(?<sub>\w+)_task-gait_acq-(?<acq>\w+)_rec-(?<rec>\d+)_nirs', 'names');
save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_nirs.mat', fparts.sub, fparts.rec)), 'data_combi')

% check events
[hdr, ~]=readoxy4(fullfile(oxy4files.folder, oxy4files.name));
stream=hdr.oxy3.Device{3}.Attributes.desc;
streamnumber=regexp(stream, '0(\d)', 'tokens');
lslstream=sprintf('lsldert%s_events', streamnumber{1}{1}); 
lsl_nirsevents=lsl_events(contains(lsl_events.type, lslstream)&contains(lsl_events.value, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}),:);
oxy4_nirsevents=nirsevents(contains({nirsevents.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
if height(lsl_nirsevents)~= length(oxy4_nirsevents)
  warning('number of events in nirs data was not as expected')
end
if any(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)>0.04)
  warning('some events of .oxy4 file differed with more than 40 msec with the lsl stream (max difference = %f)', max(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)))
else
  fprintf('maximum delay between .oxy4 file and lsl stream was %.03f seconds \n',  max(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)))
end
figure; plot([oxy4_nirsevents.sample], ones(size(oxy4_nirsevents)), '+'); title('nirs events');
% convert events to table
fsample=data_combi.fsample;
nirs_events=struct2table(oxy4_nirsevents);
onset=([oxy4_nirsevents(:).sample]'-1)/fsample;
sample=[oxy4_nirsevents(:).sample]';
onset_corrected=onset(:)-lsl_nirsevents.correction(:);
duration=zeros(length(onset),1);
type=repmat({'sync-event'}, size(onset));
value={oxy4_nirsevents(:).value}';
nirs_events=table(onset, sample, onset_corrected, duration, type, value);
% save events
save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_events.mat', fparts.sub, fparts.rec)), 'nirs_events')

% split in runs
run=[];
start_runs=oxy4_nirsevents(contains({oxy4_nirsevents.value}, 'start_run'));
stop_runs=oxy4_nirsevents(contains({oxy4_nirsevents.value}, 'stop_run'));
for r=runs
  % data
  cfg=[];
  cfg.trl=[start_runs(r).sample stop_runs(r).sample 0];
  run(r).data_nirs.data_raw=ft_redefinetrial(cfg, data_combi);
end
  
% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')

%% load motion data
% select the corresponding lsl events
fprintf('\n.........loading motion data.........\n')
lslstream='lsldert4_events';
start_runs=lsl_events(contains(lsl_events.type, lslstream) & contains(lsl_events.value, 'TTL_on'),:);
stop_runs=lsl_events(contains(lsl_events.type, lslstream) & contains(lsl_events.value, 'TTL_off'),:);

for r=runs
  fprintf('\n LOADING RUN %d \n', r)
  mvnxfile=dir(fullfile(root_dir, 'source_standard', sub, 'motion', sprintf('*_run-%02d*.mvnx', r)));
  if isempty(mvnxfile)
    error('could not find the correct motion file');
  end
  % load data
  cfg=[];
  cfg.datafile=fullfile(mvnxfile.folder, mvnxfile.name);
  data_raw=ft_preprocessing(cfg);
  % check duration of recording
  if abs((stop_runs.onset(r)-start_runs.onset(r))-data_raw.time{1}(end))>0.1
    if (stop_runs.onset(r)-start_runs.onset(r))-data_raw.time{1}(end)<0
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %d seconds longer than expected.', r, abs((stop_runs.onset(r)-start_runs.onset(r))-data_raw.time{1}(end)))
    else
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %d seconds shorter than expected.', r, abs((stop_runs.onset(r)-start_runs.onset(r))-data_raw.time{1}(end)))
    end
  else
    fprintf('maximum drift between the .mvnx file run %d and lsl stream of %.03f seconds \n', r,abs((stop_runs.onset(r)-start_runs.onset(r))-data_raw.time{1}(end)))
  end
  % check data by plotting trajectory of COM
  idx_COM=match_str(data_raw.label, {'seg_COM_centerOfMass_X','seg_COM_centerOfMass_Y','seg_COM_centerOfMass_Z'});
  figure; plot3(data_raw.trial{1}(idx_COM(1),:), data_raw.trial{1}(idx_COM(2),:),data_raw.trial{1}(idx_COM(3),:), '.'); view(2);title(sprintf('gait trajectory of run %d', r)); xlabel('X'); ylabel('Y')
  % make events
  onset = [0 data_raw.time{1}(end)]';
  onset_corrected = [0-start_runs.correction(r) data_raw.time{1}(end)-stop_runs.correction(r)]';
  duration = [0 0]';
  type={'sync-event', 'sync-event'}';
  value = {'start_run', 'stop_run'}';
  run_events=table(onset, onset_corrected, duration, type, value);
  % save data and events
  run(r).data_motion.data_raw=data_raw;
  [~, n, ~]=fileparts(mvnxfile.name);
  fparts=regexp(mvnxfile.name, 'sub-(?<sub>\w+)_task-gait_run-(?<run>\d+)_motion', 'names');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_motion.mat', fparts.sub, fparts.run)), 'data_raw');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_events.mat', fparts.sub, fparts.run)), 'run_events');
end

% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')

%% load video data
fprintf('\n.........loading video data.........\n')
acq={'mobile', 'begin', 'end'};
lslstream={'lsldert1_events', 'lsldert2_events', 'lsldert3_events'};

for a=1:length(acq)
  fprintf('LOADING CAM %s \n', acq{a})
  mp4file=dir(fullfile(root_dir, 'source_private', sub, 'video', sprintf('*acq-%s_*.mp4', acq{a})));
  if isempty(mp4file)
    error('could not find the correct video file');
  elseif length(mp4file)>1
    error('multiple recordings detected for cam %s', acq{a})
  end
  [dat, fsample, onset, value] = audio_trigger_detect(fullfile(mp4file.folder,mp4file.name));
  % check with the lslstream
  lsl_videoevents=lsl_events(contains(lsl_events.type, lslstream{a}) & contains(lsl_events.value, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}),:);
  if height(lsl_videoevents)~=length(onset)
    warning('cam %s did not contain the same amount of events as the lsl stream', acq{a})
  end
  if any(abs(onset-onset(1)-([lsl_videoevents.onset(:)]-lsl_videoevents.onset(1)))>0.1)
    warning('some events of video recording of cam %s differed with more than 100 msec with the lsl stream (max difference = %f)', max(abs(onset-onset(1)- ([lsl_videoevents.onset(:)]-lsl_videoevents.onset(1)))))
  else
    fprintf('maximum drift between video file cam %s and lsl stream of %.03f seconds \n', acq{a},max(abs(onset-onset(1)- ([lsl_videoevents.onset(:)]-lsl_videoevents.onset(1))))) 
  end
  if ~all(strcmp(value,lsl_videoevents.value))
    warning('values of the video events were different than the lsl stream')
  end
  onset_corrected=onset-lsl_videoevents.correction(:);
  duration = zeros(size(onset));
  type=repmat({'sync-event'}, size(onset));
  video_events = table(onset, onset_corrected, duration, type, value);
  % save events
  save(fullfile(root_dir, 'processed', sub, 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', sub, acq{a})), 'video_events');
end

%% save events
% sync-events
lslstream='lsldert4_events'; 
lsl_syncevents=lsl_events(contains(lsl_events.type, lslstream)&contains(lsl_events.value, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}),:);
start_runs=lsl_syncevents(contains(lsl_syncevents.value, 'start_run'),:);
stop_runs=lsl_syncevents(contains(lsl_syncevents.value, 'stop_run'),:);
% big table
events=lsl_syncevents;
events.onset=events.onset(:)-start_runs.onset(1); % relative to first start_run
events.type(:)=repmat({'lslevent'},height(events),1);
% split into runs
for r=runs
  run_events=lsl_syncevents((lsl_syncevents.onset(:)>=start_runs.onset(r) & lsl_syncevents.onset(:)<=stop_runs.onset(r)),:);
  run_events.onset=run_events.onset(:)-start_runs.onset(r);
  run_events.type(:)=repmat({'lslevent'},height(run_events),1);
  run(r).events = run_events; % store in runs
end

% FOG events
lsldert1_timestamps=lsl_correct_lsl_timestamps(REC.lsldert1_events);
lsl_FOGevents_idx=find(endsWith(REC.lsldert1_events.Data, 'FOG'));
if ~isempty(lsl_FOGevents_idx)
  onset=lsldert1_timestamps(lsl_FOGevents_idx)'-start_runs.onset(1);
  duration=zeros(size(onset));
  type=repmat({'lslevent'},length(onset),1);
  value=repmat({'FOG'}, length(onset),1);
  FOGevents=table(onset, duration, type, value);
  % add to big table
  events=[events; FOGevents];
end
 
% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_events.mat', ID)), 'events')
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')

%% save conversion script
script=mfilename('fullpath');
script_name=mfilename;
copyfile(sprintf('%s.m', script), fullfile(root_dir, 'scripts', sub, sprintf('%s_%s.m', sub, script_name)))
diary off;
