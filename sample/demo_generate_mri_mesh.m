% -------------------------------------------------------------------------
% generate mesh from .bin volume
% -------------------------------------------------------------------------

cd ..
addpath(genpath('.'))

volume_dir=['.' filesep 'volumes'];
toolbox_dir='.';
volume_name='sample_volume.bin';
mesh_filename='sample_volume_mesh.mat';

%%

volume_structure=generate_mri_mesh(volume_dir,toolbox_dir,volume_name,mesh_filename);

