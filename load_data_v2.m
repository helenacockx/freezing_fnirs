%% initialisation
close all; clear all; diary off; clc;
ft_info off
ft_notice off
ft_warning off

%%
root_dir='C:\Users\helen\Documents\freezing_fnirs\data';
source_private='F:\freezing_fnirs\source_private';
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
addpath(fullfile(root_dir, 'scripts'))
dbstop if error

ID='HC19';
sub=['sub-' ID];
sub_info=load_sub_info(ID);
runs=sub_info.runs;
num_run=sub_info.num_run;
num_block=sub_info.num_block;
% completed subjects: PD11, HC66

load_motion=0;
load_nirs=0;
load_video=1;

%% make processed folder
subfolders={'nirs', 'motion', 'video', 'stim'};
for i=1:length(subfolders)
  mkdir(fullfile(root_dir, 'processed', sub, subfolders{i}));
end
% log file
logname=fullfile(root_dir, 'processed', sub, sprintf('%s_load_data.log', sub)); diary(logname);

%% load REC variable
fprintf('\n\n.........loading lsl data.........\n')
% exceptions
switch ID
  case 'HC06' % only use rec-03
    triggerinfo=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*rec-03_triggerinfo.mat'));
  case 'PD35' % only use rec-02
      triggerinfo=dir(fullfile(root_dir, 'source_standard', sub, 'stim', '*rec-02_triggerinfo.mat'));
  otherwise
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
if load_motion
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
  data_motion=ft_preprocessing(cfg);
  % check number of channels
  if length(data_motion.label)~=782
    warning('Not all channels are present in the .mvnx file')
  end
  % check duration of recording
  duration_diff=data_motion.time{1}(end)-(stop_runs.onset(r)-start_runs.onset(r));
  if abs(duration_diff)>0.1
    if duration_diff>0
      warning('Duration of .mvnx file run %d was %.03f seconds longer than expected.', runs(r), abs(duration_diff))
    else
      warning('Duration of .mvnx file run %d was %.03f seconds shorter than expected.', runs(r), abs(duration_diff))
    end
  else
    fprintf('duration of .mvnx file run %d and lsl stream differed with %.03f seconds (.mvnx - lsl stream)\n', runs(r),duration_diff)
  end
  % check data by plotting trajectory of COM
  idx_COM=match_str(data_motion.label, {'seg_COM_centerOfMass_X','seg_COM_centerOfMass_Y','seg_COM_centerOfMass_Z'});
  figure; plot3(data_motion.trial{1}(idx_COM(1),:), data_motion.trial{1}(idx_COM(2),:),data_motion.trial{1}(idx_COM(3),:), '.'); view(2);title(sprintf('gait trajectory of run %d', runs(r))); xlabel('X'); ylabel('Y')
  % make events
  onset = [0 data_motion.time{1}(end)]';
  onset_corrected = [0 data_motion.time{1}(end)]';
  duration = [0 0]';
  type={'sync-event', 'sync-event'}';
  value = {'start_run', 'stop_run'}';
  motion_events=table(onset, onset_corrected, duration, type, value);
  % save data and events
  [~, n, ~]=fileparts(mvnxfile.name);
  fparts=regexp(mvnxfile.name, 'sub-(?<sub>\w+)_task-gait_run-(?<run>\d+)_motion', 'names');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_motion.mat', fparts.sub, fparts.run)), 'data_motion');
  save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%s_events.mat', fparts.sub, fparts.run)), 'motion_events');
end
end

%% load nirs data
if load_nirs
fprintf('\n\n.........loading nirs data.........\n')
oxy4files=dir(fullfile(root_dir, 'source_standard', sub, 'nirs', '*acq-online*.oxy4'));
r=1;
run=[];

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
  % exceptions
  switch ID
    case 'PD61'
      fprintf('first start_run event is missing in the .oxy4 file. Creating nan event. \n')
      oxy4_nirsevents=[struct('type', 'event', 'sample', nan, 'value', 'LSL start_run', 'offset', 0, 'duration', 0), oxy4_nirsevents];
    case 'PD46'
      fprintf('first 8 events are missing in the .oxy4 file. Creating nan events \n')
      missing_events=table2struct(lsl_nirsevents(1:8,:));
      oxy4_nirsevents=[struct('type', {missing_events.type}, 'sample', nan, 'value', {missing_events.value}, 'offset', 0, 'duration', 0), oxy4_nirsevents];
    case 'PD35'
      fprintf('first 2 events in .oxy4 file were not correct. Removing from list... \n')
      oxy4_nirsevents=oxy4_nirsevents(3:end);
  end
  if length(oxy4files)>1; %  multiple nirs recordings (what if first run is not 1? e.g. HC06: but first event recorded from nirs was run-05)
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
  % exceptions
  switch ID
    case 'PD61' % calculate onset of first start_run event
      onset_corrected(1)=onset(2)-(lsl_nirsevents.onset(2)-lsl_nirsevents.onset(1))-lsl_nirsevents.correction(1);
      sample(1)=round(onset_corrected(1)*data_nirs.fsample)+1; % needed to split into runs
    case 'PD46'
      onset_corrected(1:8)=onset(9)-(lsl_nirsevents.onset(9)-lsl_nirsevents.onset(1:8))-lsl_nirsevents.correction(1:8);
      sample(1:8)=round(onset_corrected(1:8)*data_nirs.fsample)+1; % needed to split into runs
  end
  nirs_events=table(onset, sample, onset_corrected,  duration, type, value);
  
  % realign events based on xsens first run
  % first resample to 60 Hz (cfr. xsens)
  cfg=[];
  cfg.resamplefs=60;
  data_nirs_rs=ft_resampledata(cfg, data_nirs);
  % select data of first run
  start_runs=nirs_events(contains(nirs_events.value, 'start_run'),:);
  stop_runs=nirs_events(contains(nirs_events.value, 'stop_run'),:);
  cfg=[];
  if strcmp(ID, 'HC76') | strcmp(ID, 'HC35') % exception: motion data of run 1 was longer than expected so not good to sync with nirs data 
    fprintf('Using run 2 instead of run 1 for syncing nirs data to motion data \n ')
    cfg.trl=[round(start_runs.onset_corrected(r+1)*data_nirs_rs.fsample)+1 round(stop_runs.onset_corrected(r+1)*data_nirs_rs.fsample)+1 0];
  else
    cfg.trl=[round(start_runs.onset_corrected(1)*data_nirs_rs.fsample)+1 round(stop_runs.onset_corrected(1)*data_nirs_rs.fsample)+1 0];
  end
  data_run=ft_redefinetrial(cfg, data_nirs_rs);
  % load xsens data and convert
  if strcmp(ID, 'HC76')|strcmp(ID, 'HC35')
    load(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%02d_motion.mat', fparts.sub, runs(r)+1)));
  else
    load(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%02d_motion.mat', fparts.sub, runs(r))));
  end
  [x, y, z]=q2e(data_motion.trial{1}(25,:), data_motion.trial{1}(26,:), data_motion.trial{1}(27,:), data_motion.trial{1}(28,:)); % convert head orientation to euler angles
  orient=rad2deg([x;y;z]); % convert from radians to degrees
  % calculate offset between nirs data and motion data
  switch ID % exception: cross correlation did not work very well
    case {'PD61', 'PD11'}
      offset_data=33;
    case 'HC33'
      offset_data=31;
    case 'HC35'
      offset_data=27;
    otherwise
      [corr, lags]=xcorr(-orient(1,:), data_run.trial{1}(102,:),30*data_run.fsample);
      [~, idx]=max(corr);
      offset_data=lags(idx);
  end
  fprintf('offset between the data streams was %.03f seconds(nirs-motion) \n', offset_data/data_run.fsample)
  fprintf('applying offset... \n')
  % update nirs_events
  nirs_events.onset_corrected=nirs_events.onset_corrected-(offset_data/data_run.fsample);
  save(fullfile(root_dir, 'processed', sub, 'nirs', sprintf('sub-%s_task-gait_rec-%s_events.mat', fparts.sub, fparts.rec)), 'nirs_events') % save events
    
  % split all data (xsens + nirs) in runs with the same time axis
  start_runs=nirs_events(contains(nirs_events.value, 'start_run'),:);
  stop_runs=nirs_events(contains(nirs_events.value, 'stop_run'),:);
  for i=1:height(start_runs)
    fprintf('\n RUN %d \n', runs(r))
    % split data
    cfg=[];
    cfg.trl=[round(start_runs.onset_corrected(i)*data_nirs_rs.fsample)+1 round(stop_runs.onset_corrected(i)*data_nirs_rs.fsample)+1 0];
    data_run_nirs=ft_redefinetrial(cfg, data_nirs_rs);
    data_run_nirs=rmfield(data_run_nirs, 'sampleinfo'); data_run_nirs.cfg=rmfield(data_run_nirs.cfg, 'trl'); % remove sampleinfo, otherwise FT still considers the original sample numbers
    % load motion data
    load(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%02d_motion.mat', fparts.sub, runs(r))));
    [x, y, z]=q2e(data_motion.trial{1}(25,:), data_motion.trial{1}(26,:), data_motion.trial{1}(27,:), data_motion.trial{1}(28,:)); % convert head orientation to euler angles
    orient=rad2deg([x;y;z]); % convert from radians to degrees
    % calculate offset between motion data and nirs data
    switch ID
      case {'PD61', 'PD11', 'HC35'}
        [corr, lags]=xcorr(-orient(1,:), data_run_nirs.trial{1}(102,:), 30);
      case 'HC35'
        [corr, lags]=xcorr(-orient(1,:), data_run_nirs.trial{1}(102,:), 100);
      otherwise
        [corr, lags]=xcorr(-orient(1,:), data_run_nirs.trial{1}(102,:), 10*data_run_nirs.fsample);
    end
    [~, idx]=max(corr);
    offset_data=lags(idx);
    fprintf('offset between the data streams was %.03f seconds (nirs-motion) \n', offset_data/data_run_nirs.fsample)
    fprintf('applying offset... \n')
    % update motion events
    load(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%02d_events.mat', fparts.sub, runs(r))));
    motion_events.onset_corrected = motion_events.onset+offset_data/data_motion.fsample;
    % save events
    save(fullfile(root_dir, 'processed', sub, 'motion', sprintf('sub-%s_task-gait_run-%02d_events.mat', fparts.sub, runs(r))), 'motion_events');
    % redefine trials
    cfg.trl=[round(motion_events.onset_corrected(1)*data_motion.fsample) round(motion_events.onset_corrected(2)*data_motion.fsample) 0];
    data_run_motion=ft_redefinetrial(cfg, data_motion);
    data_run_motion=rmfield(data_run_motion, 'sampleinfo'); data_run_motion.cfg=rmfield(data_run_motion.cfg, 'trl'); % remove sampleinfo, otherwise FT still considers the original sample numbers
    % resample
    cfg=[];
    cfg.time=data_run_nirs.time;
    data_run_motion=ft_resampledata(cfg, data_run_motion);
    % check if it worked
    [x, y, z]=q2e(data_run_motion.trial{1}(25,:), data_run_motion.trial{1}(26,:), data_run_motion.trial{1}(27,:), data_run_motion.trial{1}(28,:)); % convert head orientation to euler angles
    orient=rad2deg([x;y;z]); % convert from radians to degrees
    figure;plot(data_run_motion.time{1}, orient(1,:)*-1);
    hold on; plot(data_run_nirs.time{1}, data_run_nirs.trial{1}(102,:));
    axis([30 90 -30 30]); legend({'heading Xsens head orientation', 'heading Brite24 (24068) orientation'}); title(sprintf('xsens vs brite of run %d', runs(r)))
    % save data in run
    % exceptions 
    if strcmp(ID, 'PD90') & r==1
      fprintf('first run is empty. Do not save data in this run \n')
      r=r+1;
      continue
    elseif strcmp(ID, 'HC42') & runs(r)==3
      fprintf('Third run is empty. Do not save data in this run \n')
      r=r+1;
      continue
    elseif strcmp(ID, 'HC57') & runs(r)==2
      fprintf('During second run xsens stopped recording. Do not save data of this run.')
      r=r+1;
      continue
    end
    run(runs(r)).data_nirs.data_raw=data_run_nirs;
    run(runs(r)).data_motion.data_raw=data_run_motion;
    run(runs(r)).info.ID=ID;
    r=r+1;
  end
  
end % loop over rec

% save data
save(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)), 'run')

end
%% load video data
if load_video
fprintf('\n\n.........loading video data.........\n')
acq={'mobile', 'begin', 'end'};
lslstream={'lsldert1_events', 'lsldert2_events', 'lsldert3_events'};

  for a=1:length(acq)
    fprintf('LOADING CAM %s \n', acq{a})
%     if exist(fullfile(root_dir, 'processed', sub, 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', sub, acq{a})));
%       continue
%     end
    % exceptions
    switch ID
      case 'PD10'
        mp4file=dir(fullfile(source_private, sub, 'video', sprintf('*acq-%s_rec-01*.mp4', acq{a}))); % only use rec-01, rec-02 was on doorways
      case 'HC34' 
        if a==1
          warning('no beeps were inserted in video acq-%s of sub-%s. Skipping this video \n',acq{a}, ID)
          continue
        else
          mp4file=dir(fullfile(source_private, sub, 'video', sprintf('*acq-%s_*.mp4', acq{a})));
        end
      otherwise
        mp4file=dir(fullfile(source_private, sub, 'video', sprintf('*acq-%s_*.mp4', acq{a})));
    end
    if isempty(mp4file)
      error('could not find the correct video file');
    elseif length(mp4file)>1
      error('multiple recordings detected for cam %s', acq{a})
    end
    [dat, fsample, onset, value] = audio_trigger_detect(fullfile(mp4file.folder,mp4file.name));
    % check with the lslstream
    lsl_videoevents=lsl_events(contains(lsl_events.type, lslstream{a}) & contains(lsl_events.value, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}),:);
    % exceptions
    switch ID
      case 'PD35'
        if strcmp(acq{a}, 'mobile') | strcmp(acq{a}, 'begin')
          fprintf('first 2 events in .mp4 file were not correct. Removing from list... \n')
          onset=onset(3:end);
          value=value(3:end);
        elseif strcmp(acq{a}, 'end')
          fprintf('first 3 events in .mp4 file were not correct. Removing from list... \n')
          onset=onset(4:end);
          value=value(4:end);
          fprintf('last event is missing. Creating nan event... \n');
          onset(end+1)=nan;
          value{end+1}='stop_run';
        end
      case 'PD76'
        if strcmp(acq{a}, 'mobile')
          fprintf('first event is missing because recording was started too late. Creating nan event... \n')
          onset=[nan; onset];
          value=['start_run'; value];
        end
    end
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
    idx=find(abs(onset-onset(1)-([lsl_videoevents.onset(:)]-lsl_videoevents.onset(1)))>0.1);
    if any(idx)
      warning('some events of video recording of cam %s differed with more than 100 msec with the lsl stream (video event - lsl event):', acq{a})
      table(idx', onset(idx), value(idx), (onset(idx)-onset(1))-(lsl_videoevents.onset(idx)-lsl_videoevents.onset(1)),'VariableNames', {'idx', 'onset', 'value', 'difference'})
    else
      fprintf('maximum drift between video file cam %s and lsl stream of %.03f seconds \n', acq{a},max(abs(onset-onset(1)- ([lsl_videoevents.onset(:)]-lsl_videoevents.onset(1)))))
    end
    idx=find(~strcmp(value,lsl_videoevents.value));
    if any(idx)
      warning('some values of the video events were different than the lsl stream')
      table(idx', onset(idx), value(idx), lsl_videoevents.value(idx), 'VariableNames', {'idx', 'onset', 'value_video', 'value_lsl'})
    end
    onset_corrected=onset-lsl_videoevents.correction(:);
    % exceptions:
    if strcmp(ID, 'PD35') & a==3
      onset_corrected(end)=onset(end-1)+(lsl_videoevents.onset(end)-lsl_videoevents.onset(end-1))-lsl_videoevents.correction(end);
    elseif strcmp(ID, 'PD76') & a==1
      onset_corrected(1)=onset(2)-(lsl_videoevents.onset(2)-lsl_videoevents.onset(1))-lsl_videoevents.correction(1);
    end
    duration = zeros(size(onset));
    type=repmat({'sync-event'}, size(onset));
    video_events = table(onset, onset_corrected, duration, type, value);
    % save events
    save(fullfile(root_dir, 'processed', sub, 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', ID, acq{a})), 'video_events');
  end
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
load(fullfile(root_dir, 'processed', sub, sprintf('sub-%s_run.mat', ID)))
for r=1:length(runs)
  run_events=lsl_syncevents((lsl_syncevents.onset(:)>=start_runs.onset(r) & lsl_syncevents.onset(:)<=stop_runs.onset(r)),:);
  run_events.onset=run_events.onset(:)-start_runs.onset(r);
  run_events.type(:)=repmat({'lslevent'},height(run_events),1);
  run(runs(r)).events = run_events; % store in runs
end

% FOG events: FIXME: maybe also in run events?
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
