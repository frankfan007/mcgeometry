function dir_struct=set_default_paths(freesurfer_directory,mcx_basepath,iso2mesh_path,mctoolbox_path)

% set default paths for Monte Carlo processing
% user should be able to open and adjust paths as needed

%%

addpath(freesurfer_directory);

% working folder where processed data and volume will be saved
working_folder=mctoolbox_path;
working_folder_unix=working_folder;
if ~exist(working_folder,'dir'), mkdir(working_folder);end

% raw data folder where DCS and MRI data is stored
raw_data_basepath=[mctoolbox_path filesep 'raw' filesep];
dir_struct.dcs_file_dir=raw_data_basepath; % directory where DCS data is stored
dir_struct.mr_dir=[mctoolbox_path filesep 'volumes' filesep]; % directory where T1 is stored

% subject-based MRI volume directory, created in working folder
dir_struct.volume_dir=[working_folder filesep 'volumes' filesep];
dir_struct.volume_dir_unix=[working_folder_unix filesep 'volumes'];
if ~exist(dir_struct.volume_dir), mkdir(dir_struct.volume_dir);end

% folder to save subject-specific Monte Carlo forward simulation 
dir_struct.mc_save_dir=[working_folder filesep 'mc'];
if ~exist(dir_struct.mc_save_dir,'dir'), mkdir(dir_struct.mc_save_dir);end

% ============= software and toolbox paths ============= %

% path where MCX software is stored
dir_struct.mcx_bin=[mcx_basepath 'bin/mcx'];

% path where Monte Carlo geometry toolbox is stored
dir_struct.dcs_mc_toolbox_unix=mctoolbox_path;
dir_struct.dcs_mc_toolbox=[dcs_mc_toolbox_serverpath dir_struct.dcs_mc_toolbox_unix];

% ==================== adding paths ==================== %

addpath(genpath(iso2mesh_path));
addpath([mcx_serverpath mcx_path '/utils']);
addpath(genpath(dir_struct.dcs_mc_toolbox));
