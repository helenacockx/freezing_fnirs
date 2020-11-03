% Phase 1 (organize the recorded files)
% - Collect all files from all devices
% - Convert files in proprietary formats to open formats (oxy4->oxy3, mvn->mvnx)
% - Scan the paper lab notes and questionnaires into a pdf file
% - Rename all files following BIDS
% - Place all files in directory organization following BIDS:
% incoming2source.m
%
% ==> the source_standard data gets archived in the DAC and the source_private get archived in the Ceph storage at CNCZ?
%
% Phase 2a (add the metadata files): source2raw.m
% - Extract the timing of the events from each individual recording file as trigger-events.tsv according to BIDS
% - Determine the relative timing of all recordings  to each other and store as scans.tsv according to BIDS
% - Determine the timing of the experimental runs/blocks/trials in each of the recordings
%
% Phase 2b (minimal pre-processing)
% - Annotate the video data (freezing of gate), add it to events.tsv
% - Determine optode positions from 3D structure sensor scans
% - Remove all files with identifying information (video, 3D scans)
%
% ==> this gets archived as the â€œraw BIDS dataâ€? and shared with collaborators
%
% Phase 3 (analyze)
% - Share the anonymous data with project collaborators
% - Analyze the data
% - Publish your results and the anonymous data
%



%% Set parameters
clear all
root_dir = 'C:\Users\helen\Documents\freezing_fnirs\data';
sub='sub-HC60';
cd(root_dir)

rec_video=true
rec_nirs=[1];
nirs_offline=true;
rec_motion=[1 1 1 1];
take_motion=[1 2 3 4]; %of all recordings
run_motion=[1 2 3 4];
rec_stim=[1];

%% create source directory structure
% source_standard
mkdir(fullfile(root_dir, 'source_standard'), sub)
subdir={'nirs', 'motion', 'stim', 'questionnaires'};
for i=1:length(subdir)
    mkdir(fullfile(root_dir, 'source_standard', sub), subdir{i})
end

% source_private
mkdir(fullfile(root_dir, 'source_private'), sub)
subdir={'video', 'structuresensor'};
for i=1:length(subdir)
    mkdir(fullfile(root_dir, 'source_private', sub), subdir{i})
end

%% video folder
if rec_video
video_incom=fullfile(root_dir, 'incoming', sub, 'video');
video_outgo= fullfile(root_dir, 'source_private', sub, 'video');
% loop through all cam folders
acq={'acq-mobile', 'acq-begin', 'acq-end'};
for i=1:3
    files=dir(fullfile(video_incom, acq{i}, '*.MP4'));
    for j=1:length(files)
        filename_n=sprintf('%s_task-gait_%s_rec-%.2d_video.mp4', sub,  acq{i}, j);
        [success, message]=movefile(fullfile(files(j).folder, files(j).name),fullfile(video_outgo, filename_n));
        if success
            fprintf('\n File successfully moved to source directory: %s', filename_n);
        end
    end
end

% all files moved to source directory?
files=dir(video_incom);
if length(files)>5 %contains directory ., .., acq-mobile, acq-end, acq-begin
    warning('\n Remaining files in incoming directory \n')
else
    fprintf('\n No remaining files in incoming directory \n');
end
end

%% nirs folder
if rec_nirs
nirs_incom=fullfile(root_dir, 'incoming', sub, 'nirs');
nirs_outgo=fullfile(root_dir, 'source_standard', sub, 'nirs');
% ext={'oxy3', 'oxy4'};
ext={'oxy4'};

% move rec-prep nirs file to source directory
filename_o=sprintf('%s_rec-prep_nirs.oxy4', sub);
filename_n=sprintf('%s_rec-prep_nirs.oxy4', sub);
file=dir(fullfile(nirs_incom, filename_o));
if length(file)==1
  [success, message]=movefile(fullfile(file.folder, file.name), fullfile(nirs_outgo, filename_n));
  if success
    fprintf('\n File successfully moved to source directory: %s', filename_n);
  end
else
  warning(sprintf('No file named %s in incoming directory', filename_o));
end

% move nirs files to source directory
for i=1:length(rec_nirs)
    for e=1:length(ext)
        filename_o=sprintf('%s_rec-%.2d_nirs.%s', sub, rec_nirs(i), ext{e});
        filename_n=sprintf('%s_task-gait_acq-online_rec-%.2d_nirs.%s', sub, rec_nirs(i), ext{e});
        file=dir(fullfile(nirs_incom, filename_o));
        if length(file)==1
            [success, message]=movefile(fullfile(file.folder, file.name), fullfile(nirs_outgo, filename_n));
            if success
                fprintf('\n File successfully moved to source directory: %s', filename_n);
            end
        else
            warning(sprintf('No file named %s in incoming directory', filename_o));
        end
    end
end

% move nirs offline files to source directory
if nirs_offline
  acq={'acq-24065', 'acq-24068'};
  for a=1:2
    files=dir(fullfile(nirs_incom, acq{a}, '*.oxy3'));
    if length(rec_nirs)~=length(files)
      error('Not the same number of offline files as expected')
    end
    for i=1:length(rec_nirs)
      %         filename_o=sprintf('%s_%s_rec-%.2d_nirs.oxy3', sub, acq{a}, rec_nirs(i));
      filename_n=sprintf('%s_task-gait_%s_rec-%.2d_nirs.oxy3', sub, acq{a}, rec_nirs(i));
      [success, message]=movefile(fullfile(files(i).folder, files(i).name), fullfile(nirs_outgo, filename_n));
      if success
        fprintf('\n File successfully moved to source directory: %s', filename_n);
      end
    end
  end
end

% move DAQ screenshots to source folder
[success, message]=movefile(fullfile(nirs_incom, 'DAQ'), fullfile(nirs_outgo, 'DAQ'));
files=dir(fullfile(nirs_outgo, 'DAQ'));
if success
    fprintf('\n File successfully moved to source directory: %s consisting of %d screenshots', 'DAQ', length(files)-2);
    if length(files)-2~=8
      warning('not all 8 screenshots were captured')
    end
end


% all files moved to source directory?
files=dir(fullfile(nirs_incom, '**/*.*'));
if length(files)>8 %contains directory . and .. (x3); and acq-24065, acq-24068
    warning(sprintf('\n Remaining files in incoming directory \n'))
else
    fprintf('\n No remaining files in incoming directory \n');
end
end

%% motion folder
if rec_motion
ext={'mvn', 'mvnx'};

% move bodydimensions to source directory
motion_incom=fullfile(root_dir, 'incoming', sub, 'motion');
motion_outgo= fullfile(root_dir, 'source_standard', sub, 'motion');
files=dir(fullfile(motion_incom, '*bodydimensions.mvna'));
if length(files)==1
    [success, message]=movefile(fullfile(files(1).folder, files(1).name),fullfile(motion_outgo, files(1).name));
    if success
        fprintf('\n File successfully moved to source directory: %s', files(1).name);
    end
else
    warning('No file named *.bodydimensions.mvna or multiple files named like that');
end

% move motion files to source directory
for i=1:length(rec_motion)
    for e=1:length(ext)
        filename_o=sprintf('%s_rec-%.2d_take-%.3d.%s',sub, rec_motion(i), take_motion(i), ext{e});
        filename_n=sprintf('%s_task-gait_run-%.2d_motion.%s', sub, run_motion(i), ext{e});
        file=dir(fullfile(motion_incom, filename_o));
        if length(file)==1
            [success, message]=movefile(fullfile(file.folder, file.name),fullfile(motion_outgo,filename_n));
            if success
                fprintf('\n File successfully moved to source directory: %s', filename_n);
            end
        else
            warning('No file named %s in incoming directory',  filename_o);
        end
    end
end

% all files moved to source directory?
files=dir(motion_incom);
if length(files)>2 %contains directory . and ..
    warning('Remaining files in incoming directory \n')
else
    fprintf('No remaining files in incoming directory \n');
end
end

%% stim
if rec_stim
stim_incom=fullfile(root_dir, 'incoming', sub, 'stim');
stim_outgo= fullfile(root_dir, 'source_standard', sub, 'stim');
ext={'triggerinfo.mat', 'triggers.log', 'triggerslabrecorder.xdf'};

% move stimulation script
m_script=dir(fullfile(stim_incom, sprintf('%s_script_freezing_fnirs*.m', sub)));
if length(m_script)==1
    [success, message]=movefile(fullfile(m_script(1).folder, m_script(1).name),fullfile(stim_outgo, m_script(1).name));
    if success
        fprintf('\n File successfully moved to source directory: %s', m_script(1).name);
    end
else
    warning('\n No file named *_script_freezing_fnirs*.m in incoming directory or multiple files named like that');
end

% move stim files of practice recording to source directory
for e=1:length(ext)-1
  filename_o=sprintf('%s_rec-practice_%s', sub, ext{e});
  filename_n=sprintf('%s_task-gait_rec-practice_%s', sub, ext{e});
  file=dir(fullfile(stim_incom,filename_o));
  if length(file)==1
    [success, message]=movefile(fullfile(file.folder, file.name), fullfile(stim_outgo, filename_n));
    if success
      fprintf('\n File successfully moved to source directory: %s', filename_n);
    end
  else
    warning('\n No file named %s in incoming directory', filename_o);
  end
end

% move stim files to source directory
for i=1:length(rec_stim)
    for e=1:length(ext)
        filename_o=sprintf('%s_rec-%.3d_%s', sub, rec_stim(i), ext{e});
        filename_n=sprintf('%s_task-gait_rec-%.2d_%s', sub, rec_stim(i), ext{e});
        file=dir(fullfile(stim_incom,filename_o));
        if length(file)==1
            [success, message]=movefile(fullfile(file.folder, file.name), fullfile(stim_outgo, filename_n));
            if success
                fprintf('\n File successfully moved to source directory: %s', filename_n);
            end
        else
            warning('\n No file named %s in incoming directory', filename_o);
        end
    end
end

% all files moved to source directory?
files=dir(stim_incom);
if length(files)>2 %contains directory . and ..
    warning('\n Remaining files in incoming directory \n')
else
    fprintf('\n No remaining files in incoming directory \n');
end
end

%% structure sensor
ss_incom=fullfile(root_dir, 'incoming', sub, 'structuresensor');
ss_outgo= fullfile(root_dir, 'source_private', sub, 'structuresensor');
ext={'jpg', 'mtl', 'obj'};

% move model files to source directory
for e=1:length(ext)
    file=dir(fullfile(ss_incom, 'Model', sprintf('Model.%s', ext{e})));
    filename_n= sprintf('%s_structuresensor.%s', sub, ext{e});
    if length(file)==1
        [success, message]=movefile(fullfile(file.folder, file.name), fullfile(ss_outgo,filename_n));
        if success
            fprintf('\n File successfully moved to source directory: %s', filename_n);
        end
    else
        warning('\n No file named %s in incoming directory',  sprintf('Model.%s', ext{e}));
    end
end

% all files moved to source directory?
files=dir(fullfile(ss_incom, 'Model'));
if length(files)>2 %contains directory . and ..
    warning('\n Remaining files in incoming directory \n')
else
    fprintf('\n No remaining files in incoming directory \n');
end

%% questionnaires
quest_incom=fullfile(root_dir, 'incoming', sub, 'questionnaires');
quest_outgo= fullfile(root_dir, 'source_standard', sub, 'questionnaires');
% if contains(sub, 'PD')
% ext={'labnotes.pdf', 'inclusion.pdf', 'UPDRSIII.pdf', 'general.pdf', 'history.pfd', 'NFOGQ.pdf', 'MOCA.pdf','TMT.pdf', 'HADS.pdf', 'anxiety.pdf', 'feedback.pdf'};
% elseif contains(sub, 'HC')
%   ext={'labnotes.pdf','inclusion.pdf', 'general.pdf', 'MOCA.pdf','TMT.pdf', 'HADS.pdf', 'anxiety.pdf', 'feedback.pdf'};
% end
ext={'labnotes.pdf', 'questionnaires.pdf'};

for e=1:length(ext)
    file=dir(fullfile(quest_incom, sprintf('%s_%s', sub, ext{e})));
    filename_n=sprintf('%s_%s', sub, ext{e});
    if length(file)==1
        [success, message]=movefile(fullfile(file.folder, file.name), fullfile(quest_outgo, filename_n));
        if success
            fprintf('\n File successfully moved to source directory: %s', filename_n);
        end
    else
        warning('\n No file named %s in incoming directory',  sprintf('%s_%s', sub, ext{e}));
    end
end

% all files moved to source directory?
files=dir(fullfile(quest_incom));
if length(files)>2 %contains directory . and ..
    warning('\n Remaining files in incoming directory \n')
else
    fprintf('\n No remaining files in incoming directory \n');
end


%% save script
mkdir(fullfile(root_dir, 'scripts', sub))
script=mfilename('fullpath');
script_name=mfilename;
copyfile(sprintf('%s.m', script), fullfile(root_dir, 'scripts', sub, sprintf('%s_%s.m', sub, script_name)))
