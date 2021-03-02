function  [heel]=find_heelstrike(data, events,  run, ID, varargin);
% This function detects the first heel strike right/left after the first
% walk command of the run
%
% Use as
%   [heel]=find_heelstrike(data, events, run, ID, varargin)
%
% INPUT:
%       data         = motion data in fieldtrip data structure
%       events       = run events containing the start_block and stop_block
%       sync events
%       run          = run number of which the data was given
%       ID           = ID of which the data was given
% Additional options can be specified in key-value pairs and can be:
%       'vis'    true or false (default = true)
%
% OUTPUT
%       heel = table with the heel strike event (onset, duration, type,
%       value)
%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

fs=data.fsample;

%% find heel contacts
heel_contacts=data.trial{1}(find(contains(data.label, 'Heel_footContacts')),:);
start_block=events((strcmp(events.value, 'start_block')),:);
first_heel=find(data.time{1}>start_block.onset(1) & sum(diff([ones(2,1) heel_contacts(:,:)],1,2))==1, 1);

%% plot
if vis
  figure; plot(diff([1 heel_contacts(1,:)]));
  hold on; plot(diff([1 heel_contacts(2,:)]));
  plot(first_heel, 1, 'ro');
  ylim([-2 2]); xlim([50*fs 70*fs])
  title(sprintf('Heel strikes with detected first heel strike (sub-%s run-%02d)', ID, run))
  legend({'heel strikes left', 'heel strikes right', 'first heel strike'});
end

%% make table
onset=first_heel/fs;
duration=0;
type={'gait_event'};
value={'first_heel'};
heel=table(onset, duration, type, value);