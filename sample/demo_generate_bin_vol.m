% -------------------------------------------------------------------------
% generate binary volume from T1.mgz file
% -------------------------------------------------------------------------
%
% input:
%   dir_struct: structure with fields
%       freesurfer_dir: directory where freesurfer can be sourced
%       mr_dir: directory that contains MRI structural scan
%       volume_dir: directory to save the final volume
%   volume_name: name of volume file to be saved
%   include_fiducial: 0 or 1 flag to include fiducial marker if there is one
%   varargin:
%       if volume contains fiducial marker, varargin{1}: name of volume file with vitamin E to be saved
%
%%

dir_struct.freesurfer_dir='';
dir_struct.mr_dir=['.' filesep 'mri'];
dir_struct.volume_dir=['.' filesep 'volumes'];
