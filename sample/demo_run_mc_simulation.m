% -------------------------------------------------------------------------
% run Monte Carlo simulation
% -------------------------------------------------------------------------

% currently only supporting MCX and not tMCimg

%%

cd ..
addpath(genpath('.'))

dir_struct.input_filename=['.' filesep 'mc' filesep 'LargeSlab_MultiLyr1mm_085_mus.inp'];
dir_struct.mcx_bin='';

mc_param.gpu_number=1;
mc_param.max_detected_photons=1e5;

%%

run_mc_simulation(dir_struct,mc_param)
