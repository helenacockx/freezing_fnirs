function [data_combi, events]=offline2online(filename_online, varargin)

%% initialisation
vis=ft_getopt(varargin, 'vis', true);

[path, name, ext]=fileparts(filename_online);
fparts=regexp(name, 'sub-(?<sub>\w+)_task-gait_acq-(?<acq>\w+)_rec-(?<rec>\d+)_nirs', 'names');
filename_24065 = fullfile(path, sprintf('sub-%s_task-gait_acq-24065_rec-%s_nirs.oxy3', fparts.sub, fparts.rec));
filename_24068 = fullfile(path, sprintf('sub-%s_task-gait_acq-24068_rec-%s_nirs.oxy3', fparts.sub, fparts.rec));

% %%  initialisation 
% subID='HC57';
% subject= sprintf('sub-%s', subID);
% rec = 1;
% root_dir='C:\Users\helen\Documents\freezing_fnirs\data\';
% 
% filename_24065 = fullfile(root_dir, 'source_standard', subject, 'nirs', sprintf('%s_task-gait_acq-24065_rec-%02d_nirs.oxy3', subject, rec));
% filename_24068 = fullfile(root_dir, 'source_standard', subject, 'nirs', sprintf('%s_task-gait_acq-24068_rec-%02d_nirs.oxy3', subject, rec));
% filename_online = fullfile(root_dir, 'source_standard', subject, 'nirs',sprintf('%s_task-gait_acq-online_rec-%02d_nirs.oxy4', subject, rec));
% filename_output = fullfile(root_dir, 'raw', subject, 'nirs', sprintf('%s_task-gait_rec-%02d_nirs.snirf', subject, rec));
% 
addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\oxysoft2matlab-v1.69\priv
% addpath C:\Users\helen\Documents\MATLAB\matlab_toolboxes\fix_artinis
% cd(fullfile(root_dir, 'source_standard', subject, 'nirs'))

%% read data
% online data
cfg=[];
cfg.dataset=filename_online;
data_online=ft_preprocessing(cfg);

% offline data 24065 (==> gather metainfo form readoxy3?)
[rawOD_24065,metaInfo_24065,ADvalues_24065] = readoxy3file(filename_24065);
metaInfo_24065.OptodeTemplateID = 131;
metaInfo_24065.DPF = [6 6]; % see online file
metaInfo_24065.Position = [30 30]; % see optodetemplates.xml (DistanceScale)
metaInfo_24065.abs_K = [1.1 1.1]; % not sure what this means, taken over from online file
metaInfo_24065.abs_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24065.use_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24065.Gradient = [5 5];% not sure what this means, taken over from online file
[oxyOD2_24065, metaInfo_24065] = arrangeoxy3optodes(rawOD_24065', metaInfo_24065);

% offline data 24068
[rawOD_24068,metaInfo_24068,ADvalues_24068] = readoxy3file(filename_24068);
metaInfo_24068.OptodeTemplateID = 130;
metaInfo_24068.DPF = [6 6]; % see online file
metaInfo_24068.Position = [30 30]; % see optodetemplates.xml (DistanceScale)
metaInfo_24068.abs_K = [1.1 1.1]; % not sure what this means, taken over from online file
metaInfo_24068.abs_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24068.use_H2O = [0 0];% not sure what this means, taken over from online file
metaInfo_24068.Gradient = [5 5];% not sure what this means, taken over from online file
[oxyOD2_24068, metaInfo_24068] = arrangeoxy3optodes(rawOD_24068', metaInfo_24068);

%% map offline channels onto online channels
nsamples=min(metaInfo_24065.nbSamples, metaInfo_24068.nbSamples);
nchan=sum(strcmp(data_online.hdr.chantype, 'nirs'));
OD_offline=nan(nchan, nsamples);
% 24065
for k=1:numel(metaInfo_24065.Sys.subtemplate) % loop over subtemplates
  for i=1:length(metaInfo_24065.Sys.subtemplate(k).RxTx)
    idx=find(contains(data_online.label, metaInfo_24065.Sys.subtemplate(k).RxTx{i}));
    OD_offline(idx,:)=oxyOD2_24065.Sys.subtemplate(k).OD(1:nsamples,2*i-1:2*i)';
  end
end
% 24068
for k=1:numel(metaInfo_24068.Sys.subtemplate) % loop over subtemplates
  for i=1:length(metaInfo_24068.Sys.subtemplate(k).RxTx) % now Rx starts from Rx9 and Tx from Tx11...
    try
      parts=regexp(metaInfo_24068.Sys.subtemplate(k).RxTx{i}, 'Rx(?<RxID>\d+)-Tx(?<TxID>\d+)(?<Txsub>\D)', 'names');
      chan=sprintf('Rx%d-Tx%d%s', str2double(parts.RxID)+8, str2double(parts.TxID)+10, parts.Txsub);
    catch
      parts=regexp(metaInfo_24068.Sys.subtemplate(k).RxTx{i}, 'Rx(?<RxID>\d+)-Tx(?<TxID>\d+)', 'names');
      chan=sprintf('Rx%d-Tx%d', str2double(parts.RxID)+8, str2double(parts.TxID)+10);
    end
    idx=find(contains(data_online.label, chan));
    OD_offline(idx,:)=oxyOD2_24068.Sys.subtemplate(k).OD(1:nsamples,2*i-1:2*i)';
  end
end
% remark: this is in fact the same order as in the data_online file...

%% find offset of offline vs online files
% online files
event_online{1}=trig_detect(data_online.trial{1}, 1, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 1 (24065)
event_online{2}=trig_detect(data_online.trial{1}, 53, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 53 (24068)
for k=1:2 % find for both devices
  for i=1:length(event_online{k})-1
    if contains(event_online{k}(i).type, 'up') && contains(event_online{k}(i+1).type, 'down') && (event_online{k}(i+1).sample-event_online{k}(i).sample>50 & event_online{k}(i+1).sample-event_online{k}(i).sample<500)
      fprintf('upgoing flank detected at sample %d, verify in the figure if this is correct \n', event_online{k}(i).sample)
      break
    end
  end
  trig_online{k}=event_online{k}(i).sample;
end
if vis
  figure; subplot(2,1,1); plot(data_online.trial{1}([1 53],1:15000)'); hold on; plot(trig_online{1}, 4, 'o'); plot(trig_online{2}, 4, 'o'); title('trigger detection online data'); legend({'24065', '24068'});
end

% offline files 1 (24065) and 2 (24068)
event_offline{1}=trig_detect(OD_offline, 1, 1, 15000, 3); % find up and down going flank in the first 5 minutes in channel 1 (24065)
event_offline{2}=trig_detect(OD_offline, 53, 1, 15000, 3);% find up and down going flank in the first 5 minutes in channel 53 (24068)
for k=1:2 % find for both devices
  for i=1:length(event_offline{k})-1
    if contains(event_offline{k}(i).type, 'up') && contains(event_offline{k}(i+1).type, 'down') && (event_offline{k}(i+1).sample-event_offline{k}(i).sample>50 & event_offline{k}(i+1).sample-event_offline{k}(i).sample<500)
      fprintf('upgoing flank detected at sample %d, verify in the figure if this is correct \n', event_offline{k}(i).sample)
      break
    end
  end
  trig_offline{k}=event_offline{k}(i).sample;
end
if vis
  subplot(2,1,2); plot(OD_offline([1 53],1:15000)'); hold on; plot(trig_offline{1}, 4, 'o'); plot(trig_offline{2}, 4, 'o'); title('trigger detection offline data'); legend({'24065', '24068'});
end

% calculate offsets
offset_24065=trig_offline{1}-trig_online{1};
offset_24068=trig_offline{2}-trig_online{2};

% 
% apply offsets to the offline data 
fprintf('Applying offset of %d samples for offline file 24065 and of %d samples for offline file 24068 \n', offset_24065, offset_24068)
nsamples=size(OD_offline,2)-max(offset_24065, offset_24068);
data=nan(length(data_online.label),nsamples);
data(1:52,:)=OD_offline(1:52, offset_24065+1:offset_24065+nsamples); % 24065
data(53:96,:)=OD_offline(53:96, offset_24068+1:offset_24068+nsamples); % 24068
data(97:100,:)=ADvalues_24065(:,offset_24065+1:offset_24065+nsamples);
data(101:104,:)=ADvalues_24068(:,offset_24068+1:offset_24068+nsamples);
% and paste into online_data
data_combi=data_online;
data_combi.sampleinfo=[1 nsamples];
data_combi.time{1}=[0:1/data_combi.fsample:(nsamples-1)/data_combi.fsample];
data_combi.trial{1}=data;
data_combi.hdr.nSamples=nsamples;

%% read in events
events=ft_read_event(filename_online, 'chanindx', -1, 'type', 'event');
if vis
figure; plot([events.sample], ones(size(events)), '+');
end
% 
% xdf_file=fullfile(root_dir, 'source_standard', subject, 'stim', sprintf('%s_task-gait_rec-%02d_triggerslabrecorder.xdf', subject, rec));
% xdfevent=ft_read_event(xdf_file);
% xdfevent=xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert02'));
% xdfevent=xdfevent(3:end); % first two are initialisation events
% lsl_timestamps=[xdfevent.timestamp]-xdfevent(1).timestamp;
% 
% figure; plot(oxy_timestamps, ones(size(event)), 'b+');
% hold on; plot(lsl_timestamps, ones(size(event)), 'r+');
% figure; plot(lsl_timestamps-oxy_timestamps)
% 
%% check data
if vis
  figure; subplot(2,1,1); plot(data_online.trial{1}(3,50:end), 'b-'); hold on; plot(data_combi.trial{1}(3,50:end), 'g.'); title(sprintf('sub-%s rec-%s: 24065 vs online', fparts.sub, fparts.rec)); legend({'offline data', 'online data'});
  subplot(2,1,2); plot(data_online.trial{1}(53,50:end), 'b-'); hold on; plot(data_combi.trial{1}(53,50:end), 'g.'); title(sprintf('sub-%s rec-%s: 24068 vs online', fparts.sub, fparts.rec)); legend({'offline data', 'online data'});

  cfg                = [];
  cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
  cfg.viewmode       = 'vertical';
  cfg.channel = 'Rx*';
  cfg.event =events;
  cfg.ploteventlabels= 'colorvalue';
  cfg.plotlabels= 'yes';
  cfg.fontsize=5;
  cfg.continuous     = 'yes';
  cfg.blocksize  = 60;
  cfg.nirsscale =30;
  ft_databrowser(cfg, data_combi);
end


% %% temporary save the data
% if ~exist(fullfile(root_dir, 'processed', subject, 'nirs'))
%   mkdir(fullfile(root_dir, 'processed', subject, 'nirs'));
% end
% cd(fullfile(root_dir, 'processed', subject, 'nirs'))
% save(sprintf('%s_task-gait_rec-%02d_data_online.mat', subject, rec), 'data_online')
% save(sprintf('%s_task-gait_rec-%02d_data_combi.mat', subject, rec), 'data_combi')
% 
% %% apply correct layout
% load(fullfile(root_dir, 'source_private', subject, 'structuresensor', sprintf('%s_opto_inw.mat', subject)))
% for i=1:length(data_combi.opto.optolabel)
%   idx=find(strcmp(opto_inw.label,data_combi.opto.optolabel{i}));
%   data_combi.opto.optopos(i,:)=opto_inw.chanpos(idx,:);
% end
% data_combi.opto.unit=opto_inw.unit;
% 
% %% rename channel labels & types
% % first remove battery levels from the data and the header (channel 101 &
% % 104
% channelsel=[1:99 101:103];
% data_combi.label=data_combi.label(channelsel);
% data_combi.trial{1}=data_combi.trial{1}(channelsel,:);
% data_combi.hdr.label=data_combi.hdr.label(channelsel);
% data_combi.hdr.nChans=data_combi.hdr.nChans-2;
% data_combi.hdr.chantype=data_combi.hdr.chantype(channelsel);
% data_combi.hdr.chanunit=data_combi.hdr.chanunit(channelsel);
% 
% % rename channel labels and types
% for i=1:length(data_combi.hdr.label)
%   if strcmp(data_combi.hdr.chantype{i}, 'nirs')
%     data_combi.hdr.chantype{i}='NIRSCW';
%     data_combi.hdr.chanunit{i}='arbitrary';
%   elseif strcmp(data_combi.hdr.chantype{i}, 'AUX')
%     data_combi.hdr.chantype{i}='ACCEL'; % probably not accel!!!
%     data_combi.hrd.chanunit{i}='FIXME';
%   end
% end
    

%% data2bids
% InstitutionName             = 'Radboud University';
% InstitutionalDepartmentName = 'Donders Center for Neuroscience';
% InstitutionAddress          = 'Heyendaalseweg 135, 6525 AJ, Nijmegen, The Netherlands';
% TaskName                    = 'gait';
% TaskDescription             = 'Subjects need to walk back and forth between two taped squares on the ground in which they need to make quick half turn. Halfway they encounter a 60 cm wide doorway.';
% Instructions                = 'Follow the voice instructions. You will be asked to walk back and between the two taped squares. You step with both feet into the squares and make a quick turn (leftwards at the other side, rightwards at this side). Now and then, you will have to stand still for 30 seconds. During this period of time, you may prepare for the next block by stepping through the doorway or turning in the square. Start walking again when asked for.';
% 
% cfg = [];
% 
% cfg.dataset_description.Name                        = 'PROMPT_freezing_fnirs_pilot';
% cfg.dataset_description.DatasetType                 = 'raw';
% cfg.dataset_description.Licence                     = 'PD'; %CHECKME
% cfg.dataset_description.Authors                     = 'H.M.Cockx, R.Oostenveld, F. Nieuwhof, I.G.M. Cameron, R.J.A. van Wezel'; %CHECKME
% cfg.dataset_description.Acknowledgements            = 'A. van Setten, R. Jobse'; %FIXME;
% cfg.dataset_description.Funding                     = 'FIXME';
% cfg.dataset_description.EthicsApprovals             = 'CMO Arnhem-Nijmegen (protocol NL70915.091.19)';
% % ...
% 
% cfg.InstitutionName             = InstitutionName;
% cfg.InstitutionalDepartmentName = InstitutionalDepartmentName;
% cfg.InstitutionAddress          = InstitutionAddress;
% cfg.Manufacturer                = 'Artinis Medical Systems';
% cfg.ManufacturersModelName      = 'Brite24 + Brite24';
% cfg.DeviceSerialNumber          = '24065 + 24068'; % CHECKME
% cfg.SoftwareVersion             = 'Oxysoft 3.2.70 x64';
% cfg.nirs.CapManufacturer             = 'Artinis Medical Systems';
% cfg.nirs.CapManufacturersModelName   = 'Headcap with print (customized holes)';
% cfg.nirs.SourceType                  = 'LED';
% 
% cfg.TaskName                    = TaskName;
% cfg.TaskDescription             = TaskDescription;
% cfg.Instructions                = Instructions;
% 
% cfg.coordsystem.NIRSCoordinateSystem                            = 'CTF';
% cfg.coordsystem.NIRSCoordinateUnits                             = data_combi.opto.unit;
% cfg.coordsystem.NIRSCoordinateSystemDescription                 = 'CTF head coordinates, orientation ALS, origin between the ears';
% cfg.coordsystem.NIRSCoordinateProcessingDescription             = 'optode positions have been moved 5mm inwards to represent the locations of skin contact';
% cfg.coordsystem.FiducialsDescription                            = 'Stickers were placed on the anatomical landmarks and subsequently localized with the help of a StructureSensor'; 
% cfg.coordsystem.FiducialsCoordinates                            = sprintf('"NAS": [%s], "LPA": [%s], "RPA": [%s]', num2str(opto_inw.chanpos(strcmp(opto_inw.label,'Nz'),:),'% .2f'), num2str(opto_inw.chanpos(strcmp(opto_inw.label,'LPA'),:),'% .2f'), num2str(opto_inw.chanpos(strcmp(opto_inw.label,'RPA'),:),'% .2f'))
% cfg.coordsystem.FiducialsCoordinateSystem                       = 'CTF';
% cfg.coordsystem.FiducialsCoordinateUnits                        = data_combi.opto.unit;
% cfg.coordsystem.FiducialsCoordinateSystemDescription            = 'CTF head coordinates, orientation ALS, origin between the ears';
% 
% cfg.bidsroot = fullfile(root_dir, 'raw');
% cfg.sub = subID;
% cfg.method = 'convert'; % or convert from .oxy4 --> .snirf?
% cfg.datatype = 'nirs';
% cfg.writejson='replace'; % or merge
% cfg.writetsv='replace'; % or merge
% 
% cfg.outputfile =filename_output;
% 
% cfg.events = event; % don't try to parse the AUX channels
% % fixme: do not select TTL events or FOG events
% 
% 
% data2bids(cfg, data_combi);

%% SUBFUNCTIONS %%
function event=trig_detect(dat, chan, beginsmp, endsmp, threshold)

% select data
dat=dat(chan, beginsmp:endsmp);
% discretize the signal
dat(dat< threshold) = 0;
dat(dat>=threshold) = 1;
% convert the trigger into an event with a value at a specific sample
event=[];
begpad = dat(1);
endpad = dat(end);
difftrace = diff([begpad dat endpad]);

for j=find(difftrace~=0)
  if difftrace(j)>0
    event(end+1).type   = 'up';        % distinguish between up and down flank
    event(end  ).sample = j + beginsmp - 1;      % assign the sample at which the trigger has gone up
    event(end  ).value  = dat(j);      % assign the trigger value just _after_ going up
  elseif difftrace(j)<0
    event(end+1).type   = 'down';      % distinguish between up and down flank
    event(end  ).sample = j + beginsmp - 1;      % assign the sample at which the trigger has gone down
    event(end  ).value  = dat(j-1);    % assign the trigger value just _before_ going down
  end
end
