% -------------------------------------------------------------------------
% generate input file for slab
% -------------------------------------------------------------------------

cd ..
addpath(genpath('.'))
load LargeSlab_MultiLyr1mm_input_parameters.mat

%% volume directory and name

% setting variables
volume_direc=[pwd filesep 'volumes']; % full path to volume directory
volume_name='LargeSlab_MultiLyr1mm.bin';
input_filename=[pwd filesep 'mc' filesep 'LargeSlab_MultiLyr1mm_085_mus.inp'];

%% setting initial parameters

mc_param.volume_name=[pwd filesep 'mc' filesep volume_name];
mc_param.det_radii=[0.7 1]; % radii of each detector for forward simulation
mc_param.mus=0.85; % either one value for each tissue layer, or single value for all tissue layers
mc_param.mua=0.01; % single value for all tissue layers
mc_param.n=1.37; 
mc_param.g=0.01; 
mc_param.num_phot_launched=1e9; 

% for Monte Carlo forward simulation - please see MCX documentation
mc_param.gpu_number=0; % put 0 if using tMCimg; otherwise, check GPU number with MCX documentation 
mc_param.max_detected_photons=3e6; % max number of detected photons 

% maximum number of saved photons from the MC history file
mc_param.max_saved_photons=1e5; % number of detected photons to save

%% setting parameters from loaded input parameter mat file

mc_param.source=source_loc;
mc_param.source_unit_vec=source_unit_vec;
mc_param.detectors=detector_loc;
mc_param.vol=vol;

%% extracting number of tissue types 

% extracting number of tissue types 
fileID=fopen([volume_direc filesep volume_name]);
file_contents=fread(fileID);
fclose(fileID);

mc_param.num_tissue_layers=max(unique(file_contents));

%% check variables are consistent

% check that scattering array and number of tissue types match
if length(mc_param.mus)==1 && mc_param.num_tissue_layers>1
    mc_param.mus=repmat(mc_param.mus,1,mc_param.num_tissue_layers);
elseif length(mc_param.mus)>1 && length(mc_param.mus)~=mc_param.num_tissue_layers
    error('Number of scattering values must match number of tissues\n');
end

%% write input file

write_input_file(input_filename,mc_param);
