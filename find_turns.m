function   [turns]=find_turns(data, rot_tra, run, ID, varargin)
% This function detects the 180 degree turns based on the pelvis
% position
%
% Use as
%   [turns]=find_turns(data, rot_tra, run, ID, varargin)
%
% INPUT:
%       data         = motion data in fieldtrip data structure (not
%       reframed yet)
%       rot_tra      = rotation and translation parameters outputted from
%       ROTATION_TRANSLATION
%       run          = run number of which the data was given
%       ID           = ID of which the data was given
% Additional options can be specified in key-value pairs and can be:
%       'vis'    true or false (default = true)
%
% OUTPUT
%       turns = table with the turn events (onset, duration, type,
%       value)
%
% dependencies: REFRAME
% inspiration: Beyea et al. (2017) and script Janne/Sabine,
% https://ieeexplore.ieee.org/document/5446357, 
% https://ieeexplore.ieee.org/document/6392610?reload=true&arnumber=6392610
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample;

%% select pelvis orientation data
% cfg=[];
% cfg.channel={'seg_RightShoulder_orientation_Q0', 'seg_RightShoulder_orientation_Q1', 'seg_RightShoulder_orientation_Q2', 'seg_RightShoulder_orientation_Q3'};
% ang_shoulder=ft_selectdata(cfg, data);
% 
% cfg=[];
% cfg.channel={'seg_LeftShoulder_orientation_Q0', 'seg_LeftShoulder_orientation_Q1', 'seg_LeftShoulder_orientation_Q2', 'seg_LeftShoulder_orientation_Q3'};
% ang_shoulder=ft_selectdata(cfg, data);

cfg=[];
cfg.channel={'seg_Pelvis_orientation_Q0','seg_Pelvis_orientation_Q1', 'seg_Pelvis_orientation_Q2', 'seg_Pelvis_orientation_Q3'};
ang_pelvis=ft_selectdata(cfg, data);

% convert orientations to euler angles
[x, y, z]=q2e(ang_pelvis.trial{1}(1,:), ang_pelvis.trial{1}(2,:),ang_pelvis.trial{1}(3,:),ang_pelvis.trial{1}(4,:)); % convert orientation to euler angles
orient=rad2deg([x;y;z]);
if vis
  figure; plot(orient(3,:)); hold on;
end

% recalibrate orientation %FIXME: or use rot_angle?
start_orient=median(orient(3,[1:60*60]));
rot_tra.rot_angle;
turn_angle=orient(3,:)-start_orient;
if vis
  hold on; plot(turn_angle);
end

% sum up turning degrees (so no change from 360 to 0 degrees or vv)
idx=find(diff(turn_angle)>300|diff(turn_angle)<-300);
for i=idx
  if turn_angle(i)>turn_angle(i+1)
    turn_angle(i+1:end)=turn_angle(i+1:end)+360;
  elseif turn_angle(i)<turn_angle(i+1)
    turn_angle(i+1:end)=turn_angle(i+1:end)-360;
  end
end
hold on; plot(turn_angle);

% FIXME: low pass filter?
turn_angle_lp=ft_preproc_lowpassfilter(turn_angle, fs, 2);
if vis
  plot(turn_angle_lp)
end

%% select and reframe position data of COM
data_COM=reframe(data, 'seg_COM_centerOfMass', 'rot+tra', rot_tra);


%% find samples where orientation crosses 90 degrees, rounded at 10 degrees
crossing=find(mod(round(turn_angle_lp, -1),180)==90);
crossing=[crossing(1) crossing(find([0 diff(crossing)]>600))];% if consecutive samples, use first sample
if vis
  plot(crossing, turn_angle_lp(crossing), 'ro');
end
% find changepoints (FIX ME: adjust this to expected number of turns)
%   crossing=findchangepts(turn_angle_lp, 'MaxNumChanges', 7);
%   plot(crossing, turn_angle_lp(crossing), 'ro');

%% find start and stop of turn
turns=table('Size', [length(crossing) 4], 'VariableTypes', {'double', 'double', 'string', 'string'}, 'VariableNames', {'onset', 'duration', 'type', 'value'});
idx=1;
if contains(ID, 'PD90')
  margin=25;
elseif contains(ID, 'PD77')
  margin=20; % 20 seconds margin
else
  margin=12;
end
for i=crossing
  % exception (recording ends before 20 extra seconds)
  if strcmp(ID, 'PD90') & run==3 & i==24176
    margin=18;
  end
  turn_window_before=[i-margin*fs:i]; % window of x seconds before turning point
  turn_window_after=[i:i+margin*fs]; % window of x seconds after turning point
  if turn_angle_lp(i-margin*fs)<turn_angle_lp(i+margin*fs) % turn to the left
    turns.value(idx)='turn_180_L';
    % find local min before turn
    min_angle=min(turn_angle_lp(turn_window_before));
    loc_min=find(islocalmin(turn_angle_lp(turn_window_before), 'MinSeparation', fs/2) & turn_angle_lp(turn_window_before)<min_angle+15);
    %       loc_min=find(turn_angle_lp(turn_window_before)<turn_angle_lp(i)-80 & islocalmin(turn_angle_lp(turn_window_before), 'MinSeparation', fs/2)); %find local minimum in 40 seconds around turn
    %     loc_min=find(abs(data_COM.trial{1}(1,turn_window_before))>abs(data_COM.trial{1}(1,i))-1 & islocalmin(turn_angle(turn_window_before))); %find local minimum in 40 seconds around turn
    start_turn=turn_window_before(loc_min(end));% take last one
    % find local max after turn
    max_angle=max(turn_angle_lp(turn_window_after));
    loc_max=find(islocalmax(turn_angle_lp(turn_window_after), 'MinSeparation', fs/2) & turn_angle_lp(turn_window_after)>max_angle-15);
    %       loc_max=find(turn_angle_lp(turn_window_after)>turn_angle_lp(i)+80 & islocalmax(turn_angle_lp(turn_window_after),'MinSeparation', fs/2)); %find local minimum in 40 seconds around turn
    %     loc_max=find(abs(data_COM.trial{1}(1,turn_window_after))>abs(data_COM.trial{1}(1,i))-1  & islocalmax(turn_angle(turn_window_after))); %find local minimum in 40 seconds around turn
    stop_turn=turn_window_after(loc_max(1));
    % plot
    if vis
      plot([start_turn stop_turn], turn_angle([start_turn stop_turn]), 'g*');
    end
  elseif turn_angle(i-margin*fs)>turn_angle(i+margin*fs) % turn to the right
    turns.value(idx)='turn_180_R';
    % find local max before turn
    max_angle=max(turn_angle_lp(turn_window_before));
    loc_max=find(islocalmax(turn_angle_lp(turn_window_before), 'MinSeparation', fs/2) & turn_angle_lp(turn_window_before)>max_angle-15);
    %       loc_max=find(turn_angle_lp(turn_window_before)>turn_angle_lp(i)+80 & islocalmax(turn_angle_lp(turn_window_before),'MinSeparation', fs/2)); %find local minimum in 40 seconds around turn
%     loc_max=find(abs(data_COM.trial{1}(1,turn_window_before))>abs(data_COM.trial{1}(1,i))-1 & islocalmax(turn_angle(turn_window_before))); %find local minimum in 40 seconds around turn
    start_turn=turn_window_before(loc_max(end));
    % find local min after turn
    min_angle=min(turn_angle_lp(turn_window_after));
    loc_min=find(islocalmin(turn_angle_lp(turn_window_after), 'MinSeparation', fs/2) & turn_angle_lp(turn_window_after)<min_angle+15);
    %       loc_min=find(turn_angle_lp(turn_window_after)<turn_angle_lp(i)-80 & islocalmin(turn_angle_lp(turn_window_after), 'MinSeparation', fs/2)); %find local minimum in 40 seconds around turn
%     loc_min=find(abs(data_COM.trial{1}(1,turn_window_after))>abs(data_COM.trial{1}(1,i))-1 & islocalmin(turn_angle(turn_window_after))); %find local minimum in 40 seconds around turn
    stop_turn=turn_window_after(loc_min(1));% take first one
    if vis
      plot([start_turn stop_turn], turn_angle([start_turn stop_turn]), 'g*');
    end
  end
  % save in table
  turns.onset(idx)=start_turn/fs;
  turns.duration(idx)=(stop_turn-start_turn)/fs;
  turns.type(idx)='gait_event';
  idx=idx+1;
end

if vis
  title(sprintf('Right pelvis orientation around z-axis with detected turn events (sub-%s run-%02d)', ID, run))
  legend({'original', 'recalibrated', 'rewinded', 'filtered', 'crossing', 'start & stop of turn'});
end
