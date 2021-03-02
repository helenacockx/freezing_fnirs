%% initialisation
close all; clear all; diary off
ft_info off
ft_notice off

%%
root_dir='C:\Users\helen\Documents\freezing_fnirs\data';
source_private='F:\freezing_fnirs\source_private';
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
addpath(fullfile(root_dir, 'scripts'))
dbstop if error

ID='PD22';
sub=['sub-' ID];
sub_info=load_sub_info(ID);
runs=sub_info.runs;
num_run=sub_info.num_run;
num_block=sub_info.num_block; 
% completed subjects: PD11, HC66

%% make processed folder
subfolders={'nirs', 'motion', 'video', 'stim'};
for i=1:length(subfolders)
  mkdir(fullfile(root_dir, 'processed', sub, subfolders{i}));
end
% log file
logname=fullfile(root_dir, 'processed', sub, sprintf('%s_load_data.log', sub)); diary(logname);

%% load REC variable
fprintf('\n\n.........loading lsl data.........\n')
if strcmp(ID, 'HC06') % only use rec-03
  triggerinfo=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*rec-03_triggerinfo.mat'));
else
  triggerinfo=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*rec-0*_triggerinfo.mat'));
end
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
if length(unique(sync))~=1
  warning('The number of sync events were not equal between the lsl streams')
end
% check delay between events and calculate corrections
lsl_events=test_delay(REC);

%% load motion data
run=[];
% select the corresponding lsl events
fprintf('\n\n.........loading motion data.........\n')
lslstream='lsldert4_events';
start_runs=lsl_events(contains(lsl_events.type, lslstream) & contains(lsl_events.value, 'TTL_on'),:);
stop_runs=lsl_events(contains(lsl_events.type, lslstream) & contains(lsl_events.value, 'TTL_off'),:);

for r=1:num_run
  fprintf('\n LOADING RUN %d \n', runs(r))
  mvnxfile=dir(fullfile(root_dir, 'source_standard', sub, 'motion', sprintf('*_run-%02d*.mvnx', runs(r))));
  if isempty(mvnxfile)
    error('could not find the correct motion file');
  end
  % load data
  cfg=[];
  cfg.datafile=fullfile(mvnxfile.folder, mvnxfile.name);
  data_raw=ft_preprocessing(cfg);
  % check duration of recording
  duration_diff=data_raw.time{1}(end)-(stop_runs.onset(r)-start_runs.onset(r));
  if abs(duration_diff)>0.1
    if duration_diff>0
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %03f milliseconds longer than expected.', runs(r), abs(duration_diff)*1000)
      warning('removing the first %.03d milliseconds from the recording. Please check with video if this is correct.', abs(duration_diff)/1000)
      corr_onset=duration_diff;
    else
      warning('duration of .mvnx file run %d differed with more than 100 msec of what was expected. .mvnx was %03f seconds shorter than expected.', runs(r), abs(duration_diff)*1000)
    end
  else
    fprintf('duration of .mvnx file run %d and lsl stream differed with %03f milliseconds (.mvnx - lsl stream)\n', runs(r),duration_diff*1000)
    corr_onset=0;
  end
  % check data by plotting trajectory of COM
  idx_COM=match_str(data_raw.label, {'seg_COM_centerOfMass_X','seg_COM_centerOfMass_Y','seg_COM_centerOfMass_Z'});
  figure; plot3(data_raw.trial{1}(idx_COM(1),:), data_raw.trial{1}(idx_COM(2),:),data_raw.trial{1}(idx_COM(3),:), '.'); view(2);title(sprintf('gait trajectory of run %d', runs(r))); xlabel('X'); ylabel('Y')
  % make events
  onset = [0 data_raw.time{1}(end)]';
  onset_corrected = [corr_onset data_raw.time{1}(end)]';
  duration = [0 0]';
  type={'sync-event', 'sync-event'}';
  value = {'start_run', 'stop_run'}';
  run_events=table(onset, onset_corrected, duration, type, value);
  % save data and events
  [~, n, ~]=fileparts(mvnxfile.name);
  fparts=regexp(mvnxfile.name, 'sub-(?<sub>\w+)_task-gait_run-(?<run>\d+)_motion', 'names');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_motion.mat', fparts.sub, fparts.run)), 'data_raw');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_events.mat', fparts.sub, fparts.run)), 'run_events');
  % realign time axis and save in run
  start_run=onset_corrected(1);
  stop_run=onset_corrected(2);
  cfg=[];
  cfg.trl=[round(start_run*data_raw.fsample)+1 round(stop_run*data_raw.fsample)+1 0];
  data_run=ft_redefinetrial(cfg, data_raw);
  data_run=rmfield(data_run, 'sampleinfo'); data_run.cfg=rmfield(data_run.cfg, 'trl'); % remove sampleinfo, otherwise FT still considers the original sample numbers
  run(runs(r)).data_motion.data_raw=data_run;
end

% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')

%% load nirs data
fprintf('\n\n.........loading nirs data.........\n')
oxy4files=dir(fullfile(root_dir, 'source_standard', sub, 'nirs', '*acq-online*.oxy4'));
r=1;

for f=1:length(oxy4files)
[data_nirs, nirsevents]=offline2online_v2(fullfile(oxy4files(f).folder, oxy4files(f).name));
% temporary save data_combi
fparts=regexp(oxy4files(f).name, 'sub-(?<sub>\w+)_task-gait_acq-(?<acq>\w+)_rec-(?<rec>\d+)_nirs', 'names');
save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_nirs.mat', fparts.sub, fparts.rec)), 'data_nirs')

% check events
[hdr, ~]=readoxy4(fullfile(oxy4files(f).folder, oxy4files(f).name));
stream=hdr.oxy3.Device{3}.Attributes.desc;
streamnumber=regexp(stream, '0(\d)', 'tokens');
lslstream=sprintf('lsldert%s_events', streamnumber{1}{1}); 
lsl_nirsevents=lsl_events(contains(lsl_events.type, lslstream)&contains(lsl_events.value, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}),:);
oxy4_nirsevents=nirsevents(contains({nirsevents.value}, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
if strcmp(ID, 'PD61')
  fprintf('first start_run event is missing in the .oxy4 file. Creating nan event.')
  oxy4_nirsevents=[struct('type', 'event', 'sample', nan, 'value', 'LSL start_run', 'offset', 0, 'duration', 0), oxy4_nirsevents];
end
if length(oxy4files)>1; %  multiple nirs recordings
  fprintf('splitting lsl events according to number of runs for each nirs recording \n')
  num_run_rec=sum(contains({oxy4_nirsevents.value}, 'start_run'));
    start_rec=find(contains(lsl_nirsevents.value, 'start_run'),r);
    stop_rec=find(contains(lsl_nirsevents.value, 'stop_run'),num_run_rec+r-1);
  lsl_nirsevents=lsl_nirsevents(start_rec(end):stop_rec(end),:);
end
if height(lsl_nirsevents)~= length(oxy4_nirsevents)
  warning('number of events in nirs data was not as expected')
end
if any(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]'-oxy4_nirsevents(1).sample)/50)>0.04)
  warning('some events of .oxy4 file differed with more than 40 msec with the lsl stream (max difference = %f)', max(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]'-oxy4_nirsevents(1).sample)/50)))
else
  fprintf('maximum delay between .oxy4 file and lsl stream was %.03f seconds \n',  max(abs([lsl_nirsevents.onset(:)]-lsl_nirsevents.onset(1)-([oxy4_nirsevents(:).sample]'-oxy4_nirsevents(1).sample)/50)))
end
figure; plot([oxy4_nirsevents.sample], ones(size(oxy4_nirsevents)), '+'); title('nirs events');
  
% convert events to table
fsample=data_nirs.fsample;
nirs_events=struct2table(oxy4_nirsevents);
onset=([oxy4_nirsevents(:).sample]'-1)/fsample;
onset_corrected=onset-lsl_nirsevents.correction;
sample=[oxy4_nirsevents(:).sample]';
duration=zeros(length(onset),1);
type=repmat({'sync-event'}, size(onset));
value={oxy4_nirsevents(:).value}';
if strcmp(ID, 'PD61') % calculate onset of first start_run event
  onset_corrected(1)=onset(2)-(lsl_nirsevents.onset(2)-lsl_nirsevents.onset(1))-lsl_nirsevents.correction(1);
  sample(1)=round(onset_corrected(1)*data_nirs.fsample)+1; % needed to split into runs
end
nirs_events=table(onset, sample, onset_corrected,  duration, type, value);
% save events
% FIXME: if desynchronisation between xsens and nirs --> correct
% timestamps?
save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_events.mat', fparts.sub, fparts.rec)), 'nirs_events')

% split in runs
start_runs=nirs_events(contains(nirs_events.value, 'start_run'),:);
stop_runs=nirs_events(contains(nirs_events.value, 'stop_run'),:);
for i=1:height(start_runs)
  % split data
  cfg=[];
  cfg.trl=[start_runs.sample(i) stop_runs.sample(i) 0];
  data_run=ft_redefinetrial(cfg, data_nirs);
  % resample to 60 Hz (cfr. xsens)
  cfg=[];
  cfg.resamplefs=60;
  data_run=ft_resampledata(cfg, data_run);
  % load motion data and convert
  data_motion=run(runs(r)).data_motion.data_raw;
  [x, y, z]=q2e(data_motion.trial{1}(25,:), data_motion.trial{1}(26,:), data_motion.trial{1}(27,:), data_motion.trial{1}(28,:)); % convert head orientation to euler angles
  orient=rad2deg([x;y;z]); % convert from radians to degrees
  % calculate offset between nirs data and motion data
  [corr, lags]=xcorr(-orient(1,:), data_run.trial{1}(102,:));
  [~, idx]=max(corr);
  offset_data=lags(idx);
  fprintf('offset between the data streams was %d samples (nirs-motion) \n', offset_data)
  fprintf('applying offset... \n')
  if offset_data~=0
    % update start_runs
    start_runs.sample=start_runs.sample-round(offset_data/data_run.fsample*data_nirs.fsample)
    stop_runs.sample=stop_runs.sample-round(offset_data/data_run.fsample*data_nirs.fsample)
    % redefine trials
    cfg.trl=[start_runs.sample(i) stop_runs.sample(i) 0];
    data_run=ft_redefinetrial(cfg, data_nirs);
    % resample to 60 Hz (cfr. xsens)
    cfg=[];
    cfg.resamplefs=60;
    data_run=ft_resampledata(cfg, data_run);
    % update nirs_events...
    % FIXME
  end
%   data_run=rmfield(data_run, 'sampleinfo'); data_run.cfg=rmfield(data_run.cfg, 'trl'); % remove sampleinfo, otherwise FT still considers the original sample numers
  % check data
  figure; plot(data_motion.time{1}, orient(1,:)*-1); 
  hold on; plot(data_run.time{1}, data_run.trial{1}(102,:));
  axis([30 90 -30 30]); legend({'heading Xsens head orientation', 'heading Brite24 (24068) orientation'}); title('xsens vs brite')

  % save data in run
  run(runs(r)).data_nirs.data_raw=data_run;
  r=r+1;
end
end % loop over rec

% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')


%% load video data
fprintf('\n\n.........loading video data.........\n')
acq={'mobile', 'begin', 'end'};
lslstream={'lsldert1_events', 'lsldert2_events', 'lsldert3_events'};

for a=1:length(acq)
  fprintf('LOADING CAM %s \n', acq{a})
  mp4file=dir(fullfile(source_private, sub, 'video', sprintf('*acq-%s_*.mp4', acq{a})));
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
    % try only selecting the events of the given runs (e.g. HC06)
    start_runs_idx=find(strcmp(value, 'start_run'));
    stop_runs_idx=find(strcmp(value, 'stop_run'));
    run_values={}; run_onsets=[];
    for r=1:length(runs)
      run_values=[run_values; value(start_runs_idx(runs(r)):stop_runs_idx(runs(r)))];
      run_onsets=[run_onsets; onset(start_runs_idx(runs(r)):stop_runs_idx(runs(r)))];
    end
    if isequal(length(run_values), length(run_onsets), height(lsl_videoevents))
      fprintf('only selecting the events of the given runs \n')
      onset=run_onsets;
      value=run_values;
    else
      error('Tried to only select the events of the given runs, but still not the same amount of events as the lsl stream')
    end
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
lsl_syncevents=removevars(lsl_syncevents, 'correction');
start_runs=lsl_syncevents(contains(lsl_syncevents.value, 'start_run'),:);
stop_runs=lsl_syncevents(contains(lsl_syncevents.value, 'stop_run'),:);
% big table
events=lsl_syncevents;
events.onset=events.onset(:)-start_runs.onset(1); % relative to first start_run
events.type(:)=repmat({'lslevent'},height(events),1);
% split into runs
for r=1:length(runs)
  run_events=lsl_syncevents((lsl_syncevents.onset(:)>=start_runs.onset(r) & lsl_syncevents.onset(:)<=stop_runs.onset(r)),:);
  run_events.onset=run_events.onset(:)-start_runs.onset(r);
  run_events.type(:)=repmat({'lslevent'},height(run_events),1);
  run(runs(r)).events = run_events; % store in runs
end

% FOG events
lsldert1_timestamps=lsl_correct_lsl_timestamps(REC.lsldert1_events);
lsl_FOGevents_idx=find(endsWith(REC.lsldert1_events.Data, 'FOG'));
if length(lsl_FOGevents_idx)>1
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
ft_info on
ft_notice on
