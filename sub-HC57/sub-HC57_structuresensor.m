% this function can be used to obtain the 3D optode positions of each
% subject

%%  initialisation 
clear all;
subject= 'sub-HC57';
root_dir='C:\Users\helen\Documents\freezing_fnirs\data\';
obj_file=fullfile(root_dir, 'source_private', subject, 'structuresensor',  [subject '_structuresensor.obj']);

cd(fullfile(root_dir, 'source_private', subject, 'structuresensor'))

%% %% Load the 3D-model
head_surface = ft_read_headshape(obj_file);

% Convert the units to mm
head_surface = ft_convert_units(head_surface, 'mm');

% Visualize the mesh surface
figure; ft_plot_mesh(head_surface)

%% Specify the fiducials + optodes
% Note that I also call Cz a fiducial, which is in fact not a fiducial, but
% which we will use in a later stage to for coregistration. (If you did not
% put a sticker on the inion, you can leave Iz open)
cfg = [];
cfg.channel={};
cfg.channel{1} = 'Nz';
cfg.channel{2} = 'LPA';
cfg.channel{3} = 'RPA';
cfg.channel{4}= 'Iz';
cfg.channel{5} = 'Cz';
% change this if you have more or less receivers or transmitters (do not
% select the low output transmitters) (in this example 12 transmitters and
% 12 receivers)
cfg.channel(6:53)={'Tx1b', 'Rx2', 'Tx4', 'Rx4', 'Tx1d', 'Tx5', 'Tx3', 'Rx3', 'Tx1c', 'Rx1', 'Tx1a', 'Tx2',...
  'Tx6b', 'Rx6', 'Tx9', 'Rx8', 'Tx6d', 'Tx10', 'Tx8', 'Rx7', 'Tx6c', 'Rx5', 'Tx6a', 'Tx7',...
  'Rx9', 'Tx11a', 'Tx12', 'Rx11', 'Tx11c', 'Tx13', 'Tx15', 'Rx12', 'Tx11d', 'Tx14', 'Rx10', 'Tx11b',...
  'Rx13', 'Tx16a', 'Tx18', 'Rx15', 'Tx16c', 'Tx17', 'Tx20', 'Rx16', 'Tx16d', 'Tx19', 'Rx14', 'Tx16b'};
cfg.method = 'headshape';
opto = ft_electrodeplacement(cfg, head_surface);
% Do you want to change the anatomical labels for the axes [Y, n]? --> n
% Use "Rotate 3D" to rotate the 3D model.
% Click/unclick "Colors" to toggle the colors on and off: best is to use the color
% view for the fiducials, but the structure view for the optodes
% Use the mouse to click on fiducials/optodes and subsequently on the corresponding
% label to assign the markers. (Make sure you're not in "Rotate 3D" mode anymore!)
% If an error was made: double click on the label to remove this marker
% If ready --> press Q

%% Allign the axes of the coordinate system with the fiducial positions (ctf coordinates)
% for the mesh
clf;
cfg = [];
cfg.method        = 'fiducial';
cfg.coordsys      = 'ctf';
cfg.fiducial.nas  = opto.elecpos(1,:); %position of Nz
cfg.fiducial.lpa  = opto.elecpos(2,:); %position of LPA
cfg.fiducial.rpa  = opto.elecpos(3,:); %position of RPA
head_surface_aligned = ft_meshrealign(cfg, head_surface);

ft_plot_axes(head_surface_aligned)
ft_plot_mesh(head_surface_aligned)

% for the optodes
fid.chanpos       = [110 0 0; 0 90 0; 0 -90 0];       % CTF coordinates of the fiducials
fid.elecpos       = [110 0 0; 0 90 0; 0 -90 0];       % just like electrode positions
fid.label         = {'Nz','LPA','RPA'};    % same labels as in elec
fid.unit          = 'mm';                  % same units as mri

cfg               = [];
cfg.method        = 'fiducial';
cfg.coordsys = 'ctf';
cfg.target     = fid;                   % see above
cfg.elec          = opto;
cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
opto_aligned      = ft_electroderealign(cfg);

% Visualize the optodes on the alligned head surface
% Notice that the colorview is not that clean as the the structure view
figure;
ft_plot_mesh(head_surface_aligned)
ft_plot_sens(opto_aligned, 'elecsize', 10, 'style', 'b')

% to have the same visualization without the colors
figure;
ft_plot_mesh(removefields(head_surface_aligned, 'color'), 'tag', 'headshape', 'facecolor', 'skin', 'material', 'dull', 'edgecolor', 'none', 'facealpha', 1);
lighting gouraud
l = lightangle(0, 90);  set(l, 'Color', [1 1 1]/2)
l = lightangle(  0, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle( 90, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(180, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(270, 0); set(l, 'Color', [1 1 1]/3)
alpha 0.9
ft_plot_sens(opto_aligned, 'elecsize', 10, 'style', 'b')

% save
save([subject '_opto_aligned.mat'], 'opto_aligned')

%% Move optode inward
cfg = [];
cfg.method     = 'moveinward';
cfg.moveinward = 5; % determine distance to skin
cfg.channel = 2:length(opto_aligned.label); % do not move the nasion inward
cfg.keepchannel = true;
cfg.elec       = opto_aligned;
opto_inw = ft_electroderealign(cfg);

% visualize
figure;
ft_plot_mesh(removefields(head_surface_aligned, 'color'), 'tag', 'headshape', 'facecolor', 'skin', 'material', 'dull', 'edgecolor', 'none', 'facealpha', 1);
lighting gouraud
l = lightangle(0, 90);  set(l, 'Color', [1 1 1]/2)
l = lightangle(  0, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle( 90, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(180, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(270, 0); set(l, 'Color', [1 1 1]/3)
alpha 0.7
ft_plot_sens(opto_inw, 'elecsize', 10, 'style', 'b');

save([subject '_opto_inw.mat'], 'opto_inw')