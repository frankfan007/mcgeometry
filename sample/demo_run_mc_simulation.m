% -------------------------------------------------------------------------
% run Monte Carlo simulation
% -------------------------------------------------------------------------

% currently supporting Monte Carlo photon transport simulators
% MCX and tMCimg

%%

cd ..
addpath(genpath('.'))

%%

% make sure that this contains .inp extension if MCX and does not if tMCimg
dir_struct.input_filename=['.' filesep 'mc' filesep 'LargeSlab_MultiLyr1mm_085_mus.inp'];

% path containing executable for either MCX or tMCimg
dir_struct.montecarlo_bin='';

% only needed for MCX, not tMCimg
mc_param.gpu_number=1;
mc_param.max_detected_photons=1e5;

% 1 for MCX, 0 for tMCimg
mcx_flag=0;

%%

run_mc_simulation(dir_struct,mc_param,mcx_flag)
