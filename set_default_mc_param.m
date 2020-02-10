function mc_param=set_default_mc_param;

% for input file generation
mc_param.det_radii=[0.7 1]; % in mm
mc_param.mus=0.85; % either one value for each tissue layer, or single value for all tissue layers
mc_param.mua=0.01; % single value for all tissue layers
mc_param.n=1.37; 
mc_param.g=0.01; 
mc_param.num_phot_launched=1e9; 

% for Monte Carlo forward simulation - please see MCX documentation
mc_param.gpu_number=1; 
mc_param.max_detected_photons=3e6; % max number of detected photons 

% maximum number of saved photons from the MC history file
mc_param.max_saved_photons=1e5; % number of detected photons to save