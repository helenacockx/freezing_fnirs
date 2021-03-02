function   [doorways]=find_doorways(data, rot_tra, run, ID, varargin)
% This function detect the doorway crossings based on the COM position data
%
% Use as
%   [doorways]=find_doorways(data, rot_tra, run, ID, varargin)
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
%       doorways = table with the doorway events (onset, duration, type,
%       value)
%
% dependencies: REFRAME
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample;

%% select COM position data and reframe
data_COM=reframe(data, 'seg_COM_centerOfMass', 'rot+tra', rot_tra);
COM.pos.X = data_COM.trial{1}(1, :);
COM.pos.Y = data_COM.trial{1}(2, :);
COM.pos.Z = data_COM.trial{1}(3, :);

%% find doorway crossings 
% zero crossing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);   % Returns Zero-Crossing Indices Of Argument Vector
crossing = zci(COM.pos.X);                

% % find roi samples where x=0, rounded at 100 mm
% roi=find(round(COM.pos.X,1)==0);
% roi=roi(find([0 diff(roi)]~=1)); % if consecutive samples, use first sample
% 
% % near the roi, find sample that is closest to zero
% crossing=nan(size(roi));
% for i=1:length(roi)
%   [~, sample]=min(abs(COM.pos.X(roi(i)-5*fs:roi(i)+5*fs)));
%   crossing(i)=roi(i)-5*fs+sample;
% end

% there is never a zero crossing at the end of a run 
idx=crossing>length(COM.pos.X)-10*fs;
if sum(idx)>0
  warning('doorway crossing detected at the end of the run. Removing from event list...')
  crossing=crossing(~idx);
end



%% visualize
if vis
  figure; plot(COM.pos.X); hold on;
  plot(crossing, zeros(1, length(crossing)), 'ro');
  grid on;
  title(sprintf('center of mass trajectory along x-axis with detected doorway events (sub-%s run-%02d)', ID, run))
  legend({'COM', 'doorway event'});
end

%% events table
onset=crossing/data.fsample;
duration=zeros(size(onset));
type=repmat('gait_event', size(onset));
value=repmat('doorway', size(onset));
doorways=table(onset, duration, type, value);
