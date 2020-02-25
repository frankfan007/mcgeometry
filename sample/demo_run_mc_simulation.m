% -------------------------------------------------------------------------
% run Monte Carlo simulation
% -------------------------------------------------------------------------

% currently only supporting MCX and not tMCimg

%%

cd ..
addpath(genpath('.'))

% make sure that this contains .inp extension if MCX and does not if tMCimg
dir_struct.input_filename=['.' filesep 'mc' filesep 'LargeSlab_MultiLyr1mm_085_mus.inp'];

% path containing executable for either MCX or tMCimg
dir_struct.montecarlo_bin='';

mc_param.gpu_number=1;
mc_param.max_detected_photons=1e5;

mcx_flag=0;
%%

run_mc_simulation(dir_struct,mc_param,mcx_flag)
