function   [startstop]=find_startstop(data, events, rot_tra, run, ID, varargin)
% This function detects the 180 degree turns based on the right shoulder
% position
%
% Use as
%   [turns]=find_turns(data, events, rot_tra, run, ID, varargin)
%
% INPUT:
%       data         = motion data in fieldtrip data structure (not
%       reframed yet)
%       events       = run events containing the start_block and stop_block
%       sync events
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
% inspiration: Beyea et al. (2017) and script Janne/Sabine
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample;
%% select pelvis segment & reframe
data_pelvis=reframe(data, {'seg_Pelvis_velocity'}, 'rot', rot_tra);

% plot
if vis
  figure; a=plot(data_pelvis.trial{1}(1,:)); hold on;
end

% low pass
cfg=[];
cfg.lpfilter='yes';
% cfg.lpfreq=0.5;
cfg.lpfreq=6;
data_pelvis_lp=ft_preprocessing(cfg, data_pelvis);
if vis
  b=plot(data_pelvis_lp.trial{1}(1,:));
end

%% find start and stop walking events
% find start_walking as first point after start_block where velocity >0.1
% m/s
start_block=events(strcmp(events.value, 'start_block'),:);
onset=[]; value={}; idx=1;
for k=1:height(start_block)
  start_walking=find(abs(data_pelvis_lp.trial{1}(1,round(start_block.onset(k)*fs):end))>0.3,1,'first');
  onset(idx)=start_block.onset(k)+start_walking/fs;
  value{idx}='start_walking';
  if vis
    c=plot(onset(idx)*fs, data_pelvis_lp.trial{1}(1,round(onset(idx)*fs)), 'go');
  end
  idx=idx+1;
end

% find stop_walking as last point before stop_block where velocity <0.1
stop_block=events(strcmp(events.value, 'stop_block'),:);
for k=1:height(stop_block)
  stop_walking=find(abs(data_pelvis_lp.trial{1}(1,1:round(stop_block.onset(k)*fs+fs)))>0.3,1, 'last'); % stop_block button was sometime pressed earlier
  onset(idx)=stop_walking/fs;
  value{idx}='stop_walking';
  if vis
    d=plot(onset(idx)*fs, data_pelvis_lp.trial{1}(1,round(onset(idx)*fs)), 'ro');
  end
  idx=idx+1;
end

% create table
onset=onset';
value=value';
duration=zeros(size(onset));
type=repmat('gait_event', size(onset));
startstop=table(onset, duration, type, value);

% vis
if vis
  title(sprintf('Pelvis velocity along the x-axis with detected start and stop events (sub-%s run-%02d)', ID, run))
  legend([a,b,c,d], {'original', 'filtered', 'start walking', 'stop walking'});
end
