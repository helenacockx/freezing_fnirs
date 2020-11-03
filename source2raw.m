

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET-UP
subID = 'HC66';
sub = ['sub-' subID];
numbrun =4;
numbblock=12;

root='C:\Users\helen\Documents\freezing_fnirs\data';
source_std = fullfile(root, 'source_standard');
source_prv = fullfile(root, 'source_private');
raw = fullfile(root, 'raw');
addpath('C:\Users\helen\Documents\MATLAB\matlab_toolboxes\my_toolbox')
addpath(genpath(fullfile(root, 'scripts\')));
filelist_std = recursivedir(fullfile(source_std, sub));
filelist_prv = recursivedir(fullfile(source_prv, sub));

%% create raw directory structure
if ~exist(fullfile(raw, sub))
  mkdir(fullfile(raw, sub));
end

subdir={'nirs', 'motion', 'stim', 'video'};
for i=1:length(subdir)
  if ~exist(fullfile(raw, sub, subdir{i}))
    mkdir(fullfile(raw, sub, subdir{i}))
  end
end

%% repeating parameters
InstitutionName             = 'Radboud University';
InstitutionalDepartmentName = 'Donders Center for Neuroscience';
InstitutionAddress          = 'Heyendaalseweg 135, 6525 AJ, Nijmegen, The Netherlands';
TaskName                    = 'gait';
TaskDescription             = 'Subjects need to walk back and forth between two taped squares on the ground in which they need to make quick half turn. Halfway they encounter a 60 cm wide doorway.';
Instructions                = 'Follow the voice instructions. You will be asked to walk back and between the two taped squares. You step with both feet into the squares and make a quick turn (leftwards at the other side, rightwards at this side). Now and then, you will have to stand still for 30 seconds. During this period of time, you may prepare for the next block by stepping through the doorway or turning in the square. Start walking again when asked for.';

%%

% oxy4file = filelist(endsWith(filelist, '.oxy4'));
% for i=1:numel(oxy4file)
%   delete(oxy4file{i});
% end

%%

% mvnfile = filelist(endsWith(filelist, '.mvn'));
% for i=1:numel(mvnfile)
%   delete(mvnfile{i});
% end

%%
oxy4file = filelist_std(endsWith(filelist_std, '.oxy4'));
for i=1:numel(oxy4file)
  if contains(oxy4file{i}, 'rec-prep')
    continue
  end
  
  [data_combi, events] = offline2online(oxy4file{i});
  
  cfg = [];
  
  cfg.dataset_description.Name                        = 'PROMPT_freezing_fnirs_pilot';
  cfg.dataset_description.DatasetType                 = 'raw';
  cfg.dataset_description.Licence                     = 'PD'; %CHECKME
  cfg.dataset_description.Authors                     = 'H.M.Cockx, R.Oostenveld, F. Nieuwhof, I.G.M. Cameron, R.J.A. van Wezel'; %CHECKME
  cfg.dataset_description.Acknowledgements            = 'A. van Setten, R. Jobse'; %FIXME;
  cfg.dataset_description.Funding                     = 'FIXME';
  cfg.dataset_description.EthicsApprovals             = 'CMO Arnhem-Nijmegen (protocol NL70915.091.19)';
  ...
    
  cfg.InstitutionName             = InstitutionName;
  cfg.InstitutionalDepartmentName = InstitutionalDepartmentName;
  cfg.InstitutionAddress          = InstitutionAddress;
  cfg.Manufacturer                = 'Artinis Medical Systems';
  cfg.ManufacturersModelName      = 'Brite24 + Brite24';
  cfg.DeviceSerialNumber          = '24065 + 24068'; % CHECKME
  cfg.SoftwareVersion             = 'Oxysoft 3.2.70 x64';
  cfg.nirs.CapManufacturer             = 'Artinis Medical Systems';
  cfg.nirs.CapManufacturersModelName   = 'Headcap with print (customized holes)';
  cfg.nirs.SourceType                  = 'LED';

  cfg.TaskName                    = TaskName;
  cfg.TaskDescription             = TaskDescription;
  cfg.Instructions                = Instructions;
  
%   cfg.coordsystem.NIRSCoordinateSystem                            = 'CTF';
%   cfg.coordsystem.NIRSCoordinateUnits                             = data_combi.opto.unit;
%   cfg.coordsystem.NIRSCoordinateSystemDescription                 = 'CTF head coordinates, orientation ALS, origin between the ears';
%   cfg.coordsystem.NIRSCoordinateProcessingDescription             = 'optode positions have been moved 5mm inwards to represent the locations of skin contact';
%   cfg.coordsystem.FiducialsDescription                            = 'Stickers were placed on the anatomical landmarks and subsequently localized with the help of a StructureSensor';
%   cfg.coordsystem.FiducialsCoordinates                            = sprintf('"NAS": [%s], "LPA": [%s], "RPA": [%s]', num2str(opto_inw.chanpos(strcmp(opto_inw.label,'Nz'),:),'% .2f'), num2str(opto_inw.chanpos(strcmp(opto_inw.label,'LPA'),:),'% .2f'), num2str(opto_inw.chanpos(strcmp(opto_inw.label,'RPA'),:),'% .2f'))
%   cfg.coordsystem.FiducialsCoordinateSystem                       = 'CTF';
%   cfg.coordsystem.FiducialsCoordinateUnits                        = data_combi.opto.unit;
%   cfg.coordsystem.FiducialsCoordinateSystemDescription            = 'CTF head coordinates, orientation ALS, origin between the ears';
  
  cfg.bidsroot = raw;
  cfg.sub = subID;
  cfg.method = 'convert'; 
  cfg.datatype = 'nirs';
  cfg.writejson='replace'; % or merge
  cfg.writetsv='replace'; % or merge
  
%   [p, f, x] = fileparts(oxy3file{i});
%   cfg.outputfile = fullfile(raw, sub, 'nirs', [f x]);
%   
  cfg.events = events; 
  % fixme: do not select TTL events or FOG events
  
  data2bids(cfg, data_combi);
end

%%
mvnxfile = filelist_std(endsWith(filelist_std, '.mvnx'));
for i=1:numel(mvnxfile)
  
  cfg = [];
  
  cfg.InstitutionName             = InstitutionName;
  cfg.InstitutionalDepartmentName = InstitutionalDepartmentName;
  cfg.InstitutionAddress          = InstitutionAddress;
  cfg.Manufacturer                = 'Xsens';
  cfg.ManufacturersModelName      = 'MVN Awinda';
  cfg.DeviceSerialNumber          = '';%CHECKME
  cfg.SoftwareVersion             = 'Xsens MVN 2019.2'; %CHECKME
  
  cfg.TaskName                    = TaskName;
  cfg.TaskDescription             = TaskDescription;
  cfg.Instructions                = Instructions;

  
  cfg.bidsroot = raw;
  cfg.sub = subID;
  cfg.method = 'convert';
  cfg.datatype = 'motion';
  cfg.dataset = mvnxfile{i};
  cfg.writejson='replace';
  cfg.writetsv='replace';
  [p, f, x] = fileparts(mvnxfile{i}); 
  cfg.outputfile = fullfile(raw, sub, 'motion', f);
  
  cfg.coordsystem.MotionCoordinateSystem = 'FLU';
  cfg.coordsystem.MotionRotationRule = 'n/a'; % because using quaternions?
  cfg.coordsystem.MotionRotationOrder = 'ZXY; XZY for jointAngle(Ergo)ZXY channels';
  
  hdr = ft_read_header(cfg.dataset);
  onset = [0 (hdr.nSamples*hdr.nTrials-1)/hdr.Fs]';
  duration = [0 0]';
  value = {'start_run', 'stop_run'}';
  cfg.events = table(onset, duration, value);
  
  data2bids(cfg);
end

%%
mp4file = filelist_prv(endsWith(filelist_prv, '.mp4'));
for i=1:numel(mp4file)
  
  cfg = [];
  
  cfg.InstitutionName             = InstitutionName;
  cfg.InstitutionalDepartmentName = InstitutionalDepartmentName;
  cfg.InstitutionAddress          = InstitutionAddress;
  
  cfg.TaskName                    = TaskName;
  cfg.TaskDescription             = TaskDescription;
  cfg.Instructions                = Instructions;
  
  cfg.bidsroot = raw;
  cfg.sub = sub;
  cfg.method = 'copy';
  cfg.datatype = 'video';
  cfg.dataset = mp4file{i};
  [p, f, x] = fileparts(mp4file{i});
  cfg.outputfile = fullfile(raw, sub, 'video', [f x]);
  cfg.writejson='replace';
  cfg.writetsv='replace';
  [dat, fsample, onset, trial_type] = audio_trigger_detect(cfg.dataset, 'correction', 'none', 'visualize', true);
   
  % if no events were detected go to next file
  if isempty(onset)
      continue
  end
  duration = zeros(size(onset));
  cfg.events = table(onset, duration, trial_type);
  
  data2bids(cfg);
  
end

%% events.tsv for lsl events
xdffile = filelist_std(endsWith(filelist_std, '.xdf'));

% read in xdf events
xdfevent=ft_read_event(xdffile{1});

% checks
if length(find(contains({xdfevent.value},'TTL_aan')))~= 5*numbrun
    warning('%d start run events detected in xdf-file. This is not expected.', length(find(contains({xdfevent.value},'TTL_aan'))))
end
if length(find(contains({xdfevent.value},'TTL_uit')))~= 5*numbrun
    warning('%d stop run events detected in xdf-file. This is not expected.', length(find(contains({xdfevent.value},'TTL_uit'))))
end
if length(find(contains({xdfevent.value},'start_block')))~= 5*numbblock
    warning('%d start block events detected in xdf-file. This is not expected.', length(find(contains({xdfevent.value},'start_block'))))
end
if length(find(contains({xdfevent.value},'stop_block')))~= 5*numbblock
    warning('%d stop block events detected in xdf-file. This is not expected.', length(find(contains({xdfevent.value},'stop_block'))))
end

% find first lsl event
i=min(find(strcmp({xdfevent.value}, 'D 1 TTL_aan'))); % start of first run
% [m, i] = min([lslevent.timestamp]);

% create events.tsv file
onset=[xdfevent.timestamp]' - xdfevent(i).timestamp;
duration=zeros(size(onset));
value={xdfevent.value}';
type={xdfevent.type}';
lsl_events=table(onset, duration, type, value, 'VariableNames', {'onset', 'duration','type', 'value'});
[p, f, x] = fileparts(xdffile{1});
writetable(lsl_events, fullfile(raw, sub, 'stim', [f '_events.tsv']), 'FileType', 'text', 'Delimiter', '\t');
copyfile(xdffile{1}, fullfile(raw, sub, 'stim'));


%% copy labnotes
labnotes=filelist_std(endsWith(filelist_std, 'labnotes.pdf'));
copyfile(labnotes{1}, fullfile(raw, sub, 'labnotes'));

%% Create scans.tsv + visualize
% if multiple recordings: see organize_step2a.m in 200305_pilot/scripts

% scansfile = fullfile(raw, ['sub-' sub], ['sub-' sub '_scans.tsv']);
% scans = readtable(scansfile, 'FileType', 'text', 'Delimiter', 'tab');

% initialization
eventfile=dir('**/*_events.tsv');
filename=cell(length(eventfile),1);
acq_time = NaT(length(filename),1);
acq_time.Format='yyyy-MM-dd''T''HH:mm:ss.SSS';
cd(fullfile(raw, sub));
clear scans
j=1;

% rename all _events.tsv files to sync-events.tsv
% FIXME

% reference timeline
load(fullfile(source_std, sub, 'stim', sprintf('%s_rec-01_triggerinfo.mat', sub))) 
StartRec=REC.startRun(1).time; % start of first run as sent by matlab script
% calculate correction between sending of trigger and sending lsl event by
% lsldert0
lsldert0_events=lsl_correct_lsl_timestamps(REC.lsldert0_events);
idx=min(find(contains(REC.lsldert0_events.Data, 'TTL_aan'))); % idx of first run
corr=lsldert0_events(idx)-REC.trig_events(1).timestamp;
StartRec=StartRec+seconds(corr);
StartRec.Year=1920; StartRec.Month=01; StartRec.Day=01;
StartRec.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';

% test delay
test_delay(REC)

% labrecorder
xdffile=dir('**/*.xdf');
pieces=strsplit(xdffile(1).name, {'_', '.'});
lsl_events=readtable(fullfile(xdffile(1).folder, [sprintf('%s_', pieces{1:end-1}) 'events.tsv']), 'FileType', 'text', 'Delimiter', 'tab');
filename{j}=['stim\' xdffile(1).name];
acq_time(j)=StartRec;
j=j+1;
figure;

% nirs
oxy3file=dir('**/*.oxy3');
if numel(oxy3file)>1
  warning('multiple recordings: see organize_step2a.m in 200305_pilot/scripts')
end
lsl=1; %?
% collect all nirs events
pieces = strsplit(oxy3file(1).name, '_');
nirs_events=readtable(fullfile(oxy3file(1).folder, [sprintf('%s_', pieces{1:end-1}) 'events.tsv']) , 'FileType', 'text', 'Delimiter', 'tab');
first_run=nirs_events.onset(min(find(strcmp(nirs_events.value, 'LSL TTL_aan'))));
% collect all corresponding lsl events
lsl_nirs=lsl_events(find(strcmp(lsl_events.type,sprintf('Digital Triggers @ lsldert0%d', lsl))),:);
% internal checks
internal_checks(nirs_events, lsl_nirs)
% calculate acquision time
corr=lsl_nirs.onset(min(find(strcmp(lsl_nirs.value, 'TTL_aan')))); % correction time between lsldert00 and lsldert01
acq_time(j)=StartRec + seconds(corr) - seconds(first_run);
filename{j}=['nirs\' oxy3file(1).name];
% plot
json=jsondecodefile(fullfile(oxy3file(1).folder, [sprintf('%s_', pieces{1:end-1}) 'nirs.json']));
x_nirs=[acq_time(j); acq_time(j)+seconds(nirs_events.onset); acq_time(j)+seconds(json.RecordingDuration)];
hold on;line(x_nirs, 5*ones(size(x_nirs)), 'Marker', '+', 'Color', 'b')
j=j+1;

%video
acq={'acq-mobile', 'acq-begin', 'acq-end'};
lsl=[1 2 3];
for a=1:length(acq)
    mp4file=dir(sprintf('**/*_%s_video.mp4', acq{a}));
    if numel(oxy3file)>1
      warning('multiple recordings: see organize_step2a.m in 200305_pilot/scripts')
    end
    % collect video events
    pieces = strsplit(mp4file(1).name, '_');
    video_events=readtable(fullfile(mp4file(1).folder, [sprintf('%s_', pieces{1:end-1}) 'events.tsv']) , 'FileType', 'text', 'Delimiter', 'tab');
    first_run=video_events.onset(min(find(strcmp(video_events.value, 'start'))));
    % collect corresponding lsl events
    lsl_video=lsl_events(find(strcmp(lsl_events.type,sprintf('Digital Triggers @ lsldert0%d', lsl(a)))),:);
    % internal checks
    internal_checks(video_events, lsl_video);
    % calculate acquisition time
    corr=lsl_video.onset(min(find(strcmp(lsl_video.value, 'TTL_aan')))); % correction time between lsldert00 and other video lsl stream
    acq_time(j)=StartRec + seconds(corr) - seconds(first_run);
    filename{j}=['video\' mp4file(1).name];
    % plot
    json=jsondecodefile(fullfile(mp4file(1).folder, [sprintf('%s_', pieces{1:end-1}) 'video.json']));
    x_video=[acq_time(j); acq_time(j)+seconds(video_events.onset); acq_time(j)+seconds(json.VideoDuration)];
    hold on;line(x_video, (5-a)*ones(size(x_video)), 'Marker', '+', 'Color', 'r')
    j=j+1;
 end

% motion
mvnxfile=dir('**/*.mvnx');
lsl=4; 
% collect motion events
comb_events=[]; events=cell(length(mvnxfile),1); 
for f=1:length(mvnxfile)
    pieces = strsplit(mvnxfile(f).name, '_');
    events{f}=readtable(fullfile(mvnxfile(f).folder, [sprintf('%s_', pieces{1:end-1}) 'events.tsv']) , 'FileType', 'text', 'Delimiter', 'tab');
    events{f}.order=f*ones(size(events{f}.onset));
    comb_events=[comb_events; events{f}];
end
% collect corresponding lsl events
lsl_motion=lsl_events(find(strcmp(lsl_events.type,sprintf('Digital Triggers @ lsldert0%d', lsl)) & ismember(lsl_events.value, {'TTL_aan', 'TTL_uit'})),:);
% internal checks
internal_checks(comb_events, lsl_motion)
% calculate acquisition time
comb_events.lsl_onset=lsl_motion.onset(5:8); % run 1 & 2 are missing
for f=1:length(mvnxfile)
    events{f}=comb_events(find(comb_events.order==f), find(~contains(comb_events.Properties.VariableNames, 'order'))); %select all columns except order
    pieces = strsplit(mvnxfile(f).name, '_');
%     writetable(events{f}, fullfile(mvnxfile(f).folder, [sprintf('%s_', pieces{1:end-1}) 'events.tsv']), 'FileType', 'text', 'Delimiter', '\t');
    filename{j}=['motion\' mvnxfile(f).name];
    acq_time(j)=StartRec + seconds(events{f}.lsl_onset(1)) - seconds(events{f}.onset(1));
    x_motion=[acq_time(j); acq_time(j)+seconds(events{f}.onset(end))];
    hold on;line(x_motion, ones(size(x_motion)), 'Marker', '+', 'Color', 'g')
    j=j+1;
end
ylim([-2 8]); x_lim=xlim;
text([x_lim(1) x_lim(1) x_lim(1) x_lim(1) x_lim(1)] ,[5 4 3 2 1], {'nirs', 'video acq-mobile', 'video acq-begin', 'video acq-end', 'motion'})

% create scans.tsv
scans=table(filename, acq_time);
writetable(scans, fullfile(raw, sub, [sub '_scans.tsv']), 'FileType', 'text', 'Delimiter', '\t');


%%
%%% SUBFUNCTIONS %%%
function internal_checks(comb_events, lsl_events)
if height(comb_events)~= height(lsl_events)
    warning('%d missing events', height(lsl_events)-height(comb_events))
    type_events={'start run', 'stop run', 'start block', 'stop block', 'sync'}';
    event_numbers(1,1)=sum(contains(lsl_events.value, 'TTL_aan'))-sum(contains(comb_events.value, 'TTL_aan'));
    event_numbers(2,1)=sum(contains(lsl_events.value, 'TTL_uit'))-sum(contains(comb_events.value, 'TTL_uit'));
    event_numbers(3,1)=sum(contains(lsl_events.value, 'start_block'))-sum(contains(comb_events.value, 'start_block'));
    event_numbers(4,1)=sum(contains(lsl_events.value, 'stop_block'))-sum(contains(comb_events.value, 'stop_block'));
    event_numbers(5,1)=sum(contains(lsl_events.value, 'sync'))-sum(contains(comb_events.value, 'sync'));
    missing_events=table(type_events, event_numbers)
end
end