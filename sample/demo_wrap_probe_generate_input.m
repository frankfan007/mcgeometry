
% -------------------------------------------------------------------------
% IN PROGRESS: wrap probe around slab
% -------------------------------------------------------------------------

addpath(genpath('.'))
load LargeSlab_MultiLyr1mm_mesh.mat

%% user inputted parameters

ref_param.default_fiducial_pos=[90 75 1];
ref_param.use_default_fiducial=1; % flag: 1 to use default fiducial location and 0 to choose your own
ref_param.det_distances=[5 30]; % mm
ref_param.fiducial_pos=[]; % if not empty, will automatically use this value as fiducial position
ref_param.rotate_deg=0; % if not empty, will automatically use this value as rotational degree

%% setting input parameters from slab mesh

ref_param.vol=permuted_vol;

ref_param.all_node=all_node;
ref_param.all_elem=all_elem;
ref_param.all_face=all_face;

%% wrap probe

ref_param=wrap_probe(ref_param);

%%

save LargeSlab_MultiLyr1mm_input_parameters.mat ref_param