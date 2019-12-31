
function [dir_struct,ref_param,mc_param]=set_parameters_for_given_volume(dir_struct,ref_param,mc_param,volume_cfg)

% setting parameters for given volume
% input:
%   dir_struct: structure with fields
%       dcs_mc_toolbox: path to mc geometry toolbox
%       mc_save_dir: path to save results and subject-specific Monte Carlo forward simulation
%   ref_param: structure with fields
%       has_fiducial: 0 or 1 flag denoting whether scan has fiducial
%   mc_param: 
%       inp_filename: input filename without the extension
%       volume_name_noext: if subject-specific volume, name of volume without extension
%   volume_cfg: structure denoting which volume to fit with
%       multi_layer_slab: 0 or 1 flag
%       multi_layer_head: 0 or 1 flag
%       subj_specific_mri: 0 or 1 flag

% output:
%   dir_struct: added fields include
%       input_filename: full path to input filename with extension
%   ref_param: added fields include
%       default_fiducial_pos: default fiducial position
%   mc_param: added fields include
%       volume_name: name of volume
%       meshvol_name: name of mesh that accompanies volume

%%

if volume_cfg.multi_layer_slab + volume_cfg.multi_layer_head + volume_cfg.subj_specific_mri < 1
    error('Please specify at least one option for volume configurations\n');
elseif volume_cfg.multi_layer_slab + volume_cfg.multi_layer_head + volume_cfg.subj_specific_mri >1
    error('Please specify only one option for volume configurations\n');
end

if volume_cfg.multi_layer_slab
    ref_param.default_fiducial_pos=[90 75 1];
    mc_param.volume_name='LargeSlab_MultiLyr1mm.bin';
    mc_param.meshvol_name='LargeSlab_MultiLyr1mm_mesh.mat';
    dir_struct.input_filename=[dir_struct.dcs_mc_toolbox filesep 'mc' filesep mc_param.inp_filename '.inp']; 
elseif volume_cfg.multi_layer_head
    ref_param.default_fiducial_pos=[83.30 127.40 200.60];
    mc_param.volume_name='stairhead.bin';
    mc_param.meshvol_name='stairhead_mesh.mat';
    dir_struct.input_filename=[dir_struct.dcs_mc_toolbox filesep 'mc' filesep mc_param.inp_filename '.inp']; 
elseif volume_cfg.subj_specific_mri
    ref_param.default_fiducial_pos=[83.30 127.40 200.60];
    if ref_param.has_fiducial
        mc_param.volume_name_vite=[mc_param.volume_name_noext '_vite.bin'];
        mc_param.meshvol_name_vite=[mc_param.volume_name_noext '_vite_mesh.mat'];
    end
    mc_param.volume_name=[mc_param.volume_name_noext '.bin'];
    mc_param.meshvol_name=[mc_param.volume_name_noext '_mesh.mat']; 
    dir_struct.input_filename=[dir_struct.mc_save_dir filesep mc_param.inp_filename '.inp']; 
end

