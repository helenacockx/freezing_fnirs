%% initialisation
close all; clear all
ID='HC90';
sub=['sub-' ID];
root_dir='C:\Users\helen\Documents\freezing_fnirs\data';
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
addpath(fullfile(root_dir, 'scripts'))
dbstop if error
% completed subjects: PD11, HC66

runs=[5:8];
num_run=8;
num_block=12; 

%nirs
delay_devices=0;

% make processed folder
subfolders={'nirs', 'motion', 'video', 'stim'};
for i=1:length(subfolders)
  mkdir(fullfile(root_dir, 'processed', sub, subfolders{i}));
end
% log file
logname=fullfile(root_dir, 'processed', sub, sprintf('%s_load_data.log', sub)); diary(logname);


%% load xdf files
fprintf('\n.........loading lsl data.........\n')
xdffiles=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*triggerslabrecorder.xdf'));

% read in events
xdfevent=ft_read_event(fullfile(xdffiles.folder, xdffiles.name));

% check number of events
if ~isequal(sum(contains({xdfevent.value},'TTL_on')),sum(contains({xdfevent.value},'TTL_off')),sum(contains({xdfevent.value},'start_run')),sum(contains({xdfevent.value},'stop_run')), 5*num_run)
  warning('The number of start_runs/TTL_on and stop_runs/TTL_off were not as expected')
end
if ~isequal(sum(contains({xdfevent.value},'start_block')),sum(contains({xdfevent.value},'stop_block')), 5*num_block)
  warning('The number of start_blocks and stop_blocks were not as expected')
end
if ~isequal(sum(contains({xdfevent.value}, 'sync') & strcmp({xdfevent.type}, 'Digital Triggers @ lsldert00')),sum(contains({xdfevent.value}, 'sync') & strcmp({xdfevent.type}, 'Digital Triggers @ lsldert01')),sum(contains({xdfevent.value}, 'sync') & strcmp({xdfevent.type}, 'Digital Triggers @ lsldert02')),sum(contains({xdfevent.value}, 'sync') & strcmp({xdfevent.type}, 'Digital Triggers @ lsldert03')), sum(contains({xdfevent.value}, 'sync') & strcmp({xdfevent.type}, 'Digital Triggers @ lsldert04')))
  warning('The number of sync events were not equal between the lsl streams')
end
% check delay between events
test_xdfdelay(xdfevent);

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
streamnumber=regexp(stream, '\d+', 'match');
lslstream=sprintf('lsldert%s', streamnumber{1}); 
lsl_nirsevents=xdfevent(contains({xdfevent.type}, lslstream)&contains({xdfevent.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
oxy4_nirsevents=nirsevents(contains({nirsevents.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
if length(lsl_nirsevents)~= length(oxy4_nirsevents)
  warning('number of events in nirs data was not as expected')
end
if any(abs([lsl_nirsevents(:).timestamp]-lsl_nirsevents(1).timestamp-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)>0.04)
  warning('some events of .oxy4 file differed with more than 40 msec with the lsl stream (max difference = %f)', max(abs([lsl_nirsevents(:).timestamp]-lsl_nirsevents(1).timestamp-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)))
else
  fprintf('maximum delay between .oxy4 file and lsl stream was %.03f seconds \n',  max(abs([lsl_nirsevents(:).timestamp]-lsl_nirsevents(1).timestamp-([oxy4_nirsevents(:).sample]-oxy4_nirsevents(1).sample)/50)))
end
figure; plot([oxy4_nirsevents.sample], ones(size(oxy4_nirsevents)), '+'); title('nirs events');
% save events
save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_events.mat', fparts.sub, fparts.rec)), 'data_combi')


% split in runs
run=[];
start_runs=oxy4_nirsevents(contains({oxy4_nirsevents.value}, 'start_run'));
stop_runs=oxy4_nirsevents(contains({oxy4_nirsevents.value}, 'stop_run'));
fsample=data_combi.fsample;
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
lslstream='Digital Triggers @ lsldert04';
start_runs=ft_filter_event(xdfevent, 'type', lslstream, 'value', 'TTL_on');
stop_runs=ft_filter_event(xdfevent, 'type', lslstream, 'value', 'TTL_off');

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
  if abs((stop_runs(r).timestamp-start_runs(r).timestamp)-data_raw.time{1}(end))>0.1
    if (stop_runs(r).timestamp-start_runs(r).timestamp)-data_raw.time{1}(end)<0
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %d seconds longer than expected.', r, abs((stop_runs(r).timestamp-start_runs(r).timestamp)-data_raw.time{1}(end)))
    else
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %d seconds shorter than expected.', r, abs((stop_runs(r).timestamp-start_runs(r).timestamp)-data_raw.time{1}(end)))
    end
  else
    fprintf('maximum drift between the .mvnx file run %d and lsl stream of %.03f seconds \n', r,abs((stop_runs(r).timestamp-start_runs(r).timestamp)-data_raw.time{1}(end)))
  end
  % check data by plotting trajectory of COM
  idx_COM=match_str(data_raw.label, {'seg_COM_centerOfMass_X','seg_COM_centerOfMass_Y','seg_COM_centerOfMass_Z'});
  figure; plot3(data_raw.trial{1}(idx_COM(1),:), data_raw.trial{1}(idx_COM(2),:),data_raw.trial{1}(idx_COM(3),:), '.'); view(2);title(sprintf('gait trajectory of run %d', r)); xlabel('X'); ylabel('Y')
  % make events
  onset = [0 data_raw.time{1}(end)]';
  duration = [0 0]';
  type={'sync-event', 'sync-event'}';
  value = {'start_run', 'stop_run'}';
  run_events=table(onset, duration, type, value);
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
lslstream={'Digital Triggers @ lsldert01', 'Digital Triggers @ lsldert02', 'Digital Triggers @ lsldert03'};

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
  lsl_videoevents=xdfevent(contains({xdfevent.type}, lslstream{a}) & contains({xdfevent.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
  if length(lsl_videoevents)~=length(onset)
    warning('cam %s did not contain the same amount of events as the lsl stream', acq{a})
  end
  if any(abs(onset-onset(1)- ([lsl_videoevents(:).timestamp]'-lsl_videoevents(1).timestamp))>0.1)
    warning('some events of video recording of cam %s differed with more than 100 msec with the lsl stream (max difference = %f)', max(abs(onset-onset(1)- ([lsl_videoevents(:).timestamp]'-lsl_videoevents(1).timestamp))))
  else
    fprintf('maximum drift between video file cam %s and lsl stream of %.03f seconds \n', acq{a},max(abs(onset-onset(1)- ([lsl_videoevents(:).timestamp]'-lsl_videoevents(1).timestamp)))) 
  end
  if ~all(strcmp(value,{lsl_videoevents.value}'))
    warning('values of the video events were different than the lsl stream')
  end
  duration = zeros(size(onset));
  type=repmat({'sync-event'}, size(onset));
  video_events = table(onset, duration, type, value);
  % save events
  save(fullfile(root_dir, 'processed', sub, 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', sub, acq{a})), 'video_events');
end

%% save events
% sync-events
lslstream='lsldert04'; 
lsl_syncevents=xdfevent(contains({xdfevent.type}, lslstream)&contains({xdfevent.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
start_runs=lsl_syncevents(contains({lsl_syncevents.value}, 'start_run'));
stop_runs=lsl_syncevents(contains({lsl_syncevents.value}, 'stop_run'));
% big table
events=struct2table(lsl_syncevents);
events.onset=events.timestamp(:)-start_runs(1).timestamp; % relative to first start_run
events.type(:)=repmat({'lslevent'},height(events),1);
events.duration=zeros(height(events),1);
events=events(:,[7 6 1 2]); % only retain onset-duration-type-value
% split into runs
for r=runs
  run_events=ft_filter_event(lsl_syncevents, 'mintimestamp', start_runs(r).timestamp, 'maxtimestamp', stop_runs(r).timestamp);
  run_events=struct2table(run_events);
  run_events.onset=run_events.timestamp(:)-start_runs(r).timestamp;
  run_events.type(:)=repmat({'lslevent'},height(run_events),1);
  run_events.duration=zeros(height(run_events),1);
  run_events=run_events(:,[7 6 1 2]); % only retain onset-duration-type-value
  run(r).events = run_events; % store in runs
end

% FOG events
lslstream='lsldert01';
lsl_FOGevents=xdfevent(contains({xdfevent.type}, lslstream)&contains({xdfevent.value}, {'FOG'}));
if ~isempty(lsl_FOGevents)
  FOGevents=struct2table(lsl_FOGevents);
  FOGevents.onset=FOGevents.timestamp(:)-start_runs(1).timestamp; % relative to first start_run
  FOGevents.type(:)=repmat({'lslevent'},height(FOGevents),1);
  FOGevents.duration=zeros(height(FOGevents),1);
  FOGevents=FOGevents(:,[7 6 1 2]); % only retain onset-duration-type-value
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
