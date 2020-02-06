% -------------------------------------------------------------------------
% run Monte Carlo simulation
% -------------------------------------------------------------------------

% input:
%   dir_struct: structure with fields
%       input_filename: full filename of input file
%       mcx_bin: path containing mcx executable
%   mc_param:
%       gpu_number: number of GPU to use
%       max_detected_photons: maximum number of detected photons

dir_struct.input_filename=['.' filesep 'mc' filesep 'LargeSlab;

run_mc_simulation(dir_struct,mc_param)