% -------------------------------------------------------------------------
% IN PROGRESS: generate input file for slab
% -------------------------------------------------------------------------

addpath(genpath('.'))
load LargeSlab_MultiLyr1mm_085_mus_input_parameters.mat

%%

% setting variables
volume_direc=['.' filesep 'volumes'];
volume_name='LargeSlab_MultiLyr1mm.bin';

% extracting number of tissue types 
fileID=fopen([volume_direc filesep volume_name]);
file_contents=fread(fileID);
fclose(fileID);

mc_param.num_tissue_layers=max(unique(file_contents));

% moving some stuff over
mc_param.source=ref_param.all_locations(1,:);
mc_param.source_unit_vec=ref_param.source_unit_vec;
mc_param.detectors=ref_param.all_locations(2:end,:);
mc_param.fiducial_point=ref_param.fiducial_point;
mc_param.vol=ref_param.vol;

% check

if length(mc_param.mus)==1 && mc_param.num_tissue_layers>1
    mc_param.mus=repmat(mc_param.mus,1,mc_param.num_tissue_layers);
elseif length(mc_param.mus)>1 && length(mc_param.mus)~=mc_param.num_tissue_layers
    error('Number of scattering values must match number of tissues\n');
end

% write input file
% dir_struct.input_filename=[dir_struct.mc_save_dir filesep mc_param.input_file_basename '.inp'];
write_input_file(dir_struct.input_filename,mc_param);
