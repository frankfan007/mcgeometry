% -------------------------------------------------------------------------
% generate binary volume from T1.mgz file
% -------------------------------------------------------------------------

%%

dir_struct.freesurfer_dir='';
dir_struct.mr_dir=['.' filesep 'mri'];
dir_struct.volume_dir=['.' filesep 'volumes'];

%%

generate_bin_vol(dir_struct,'sample_volume.bin',0);
