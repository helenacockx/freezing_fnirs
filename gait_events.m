%% script to create the events to import in ELAN & to make table of offsets of videos to synchronize in ELAN
root_dir='C:\Users\helen\Documents\freezing_fnirs\data';
video_dir='\\dcn-srv.science.ru.nl\dcn\biophysics\prompt\freezing_fnirs\data';

addpath(fullfile(root_dir,'scripts'))
ft_info off
ft_notice off
ft_warning off

dbstop if error

ID={'HC76'};

%% create table with  video offsets to synchronize cams in ELAN
% create empty table or load last one
try
  video_offsets=readtable(fullfile(video_dir, 'processed', 'video_offsets.xlsx'));
catch
  % create new table
  sub=cell(198,1);
  for i=1:99
    sub{i}=sprintf('PD%02d', i);
  end
  for i=100:198
    sub{i}=sprintf('HC%02d', i-99);
  end
  offset_mobile=nan(198,1);
  offset_begin=nan(198,1);
  offset_end=nan(198,1);
  video_offsets=table(sub, offset_mobile, offset_begin, offset_end);
end

% fill in table
for j=1:length(ID)
  idx=find(strcmp(video_offsets.sub, ID(j)));
  load(fullfile(root_dir, 'processed', ['sub-' ID{j}], 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', ID{j}, 'mobile')));
  offset=video_events.onset_corrected(find(strcmp(video_events.value, 'start_run'),1));
  video_offsets.offset_mobile(idx)=round(offset*1000); % convert to msec to input in ELAN
  load(fullfile(root_dir, 'processed',['sub-' ID{j}], 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', ID{j}, 'begin')));
  offset=video_events.onset_corrected(find(strcmp(video_events.value, 'start_run'),1));
  video_offsets.offset_begin(idx)=round(offset*1000);
  load(fullfile(root_dir, 'processed', ['sub-' ID{j}], 'video', sprintf('sub-%s_task-gait_acq-%s_events.mat', ID{j}, 'end')));
  offset=video_events.onset_corrected(find(strcmp(video_events.value, 'start_run'),1));
  video_offsets.offset_end(idx)=round(offset*1000);
end

% save table as .xls file
writetable(video_offsets, fullfile(video_dir, 'processed', 'video_offsets.xlsx'));

%% collect motion events to import in ELAN
close all
for j=1:length(ID)
  % first select the gait tasks and FOG events from the event table
  load(fullfile(root_dir, 'processed', ['sub-' ID{j}], sprintf('sub-%s_events.mat', ID{j})));
  onset=events.onset(find(strcmp(events.value, 'start_run')));
  duration=events.onset(find(strcmp(events.value, 'stop_run')))-onset(:);
  type=repmat({'gait_task'}, size(onset));
  value=repmat({'FOG_Course'}, size(onset));
  gait_task=table(onset,duration, type, value);
  % add the FOG events
  FOG_events=events(strcmp(events.type, 'lslevent') & strcmp(events.value, 'FOG'),:);
  motion_events=[gait_task; FOG_events];
  
  % now load the motion data of each run and create gait events
  load(fullfile(root_dir, 'processed', ['sub-' ID{j}], sprintf('sub-%s_run.mat', ID{j})));
  for i=1:length(run)
    if isempty(run(i).data_motion)
      continue
    end
    % exceptions
    %   if strcmp(ID, 'PD10') & i==2
    %     warning('first part of the run contains a lot of drift. Not using this part of the data...')
    %     run(i).data_motion.data_raw.trial{1}(:,1:7030)=nan;
    %   end
    % rotate and translate position data in order to have the lenght of the corridor along the x-axis
    [rot_tra] = rotation_translation(run(i).data_motion.data_raw, i, ID{j});
    
    % define doorway events
    [doorways]=find_doorways(run(i).data_motion.data_raw, rot_tra, i, ID{j});
    
    % find turning
    [turns]=find_turns(run(i).data_motion.data_raw, rot_tra, i, ID{j});
    
    % find start and stop walking
    [startstop]=find_startstop(run(i).data_motion.data_raw, run(i).events, rot_tra, i, ID{j});
    
    % find first heel strike (to check synchronisation between xsens and video
    [heel]=find_heelstrike(run(i).data_motion.data_raw, run(i).events,  i, ID{j});
    
    % save gait events in one table
    gait_events=[heel; doorways; turns; startstop];
    gait_events=sortrows(gait_events); %sort ons ascending onset times  
    
    % add offsets of each run
    gait_events.onset=gait_events.onset+onset(i);
      
    % add to motion event table
    motion_events=[motion_events; gait_events];
      
  end
  
  % change order of types
  motion_events=sortrows(motion_events, {'type', 'onset'});
  
  % add a column end_time and add 1 second for all events with 0 s. duration
  motion_events.end_time=motion_events.onset+motion_events.duration;
  zeros=find(motion_events.duration==0);
  motion_events.end_time(zeros)=motion_events.end_time(zeros)+1;
  
  % exception
  if strcmp(ID{j}, 'PD76')
    % acq-mobile was started 25 seconds too late --> shift annotations to
    % be able to sync them in ELAN (! shift back when reading in
    % annotations)
    warning('acq-mobile was started 25.361 seconds too late. Shifting annotations to sync them in ELAN')
    motion_events.onset=motion_events.onset-25.361; 
    motion_events.end_time=motion_events.end_time-25.361;
  end
  
  % save motion events for video annotations
  display(motion_events)
  response=input('Do you want to save these motion events? 1/0 \n');
  if response
    fprintf('writing .tsv file... \n')
%     save(fullfile(video_dir, 'processed', 'motion_events', sprintf('sub-%s_motion-events.mat', ID{j})), 'motion_events');
    writetsv(motion_events, fullfile(video_dir, 'processed', 'motion_events', sprintf('sub-%s_motion-events.tsv', ID{j})));
  end
end
