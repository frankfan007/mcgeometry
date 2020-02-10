% -------------------------------------------------------------------------
% generate mesh from .bin volume
% -------------------------------------------------------------------------

% generates mesh from volume file
%
% input:
%   volume_dir: directory where volume is stored
%   toolbox_dir: directory where toolbox functions are stored
%   volume_name: name of volume
%   mesh_filename: name of mesh accompanying the volume
%   varargin: if volume contains fiducial
%       varargin{1}: name of volume with fiducial markers included
%       varargin{2}: name of mesh accompanying volume with fiducial markers included
%
% output:
%   volume_structure: structure with fields
%       vol: volume, dimension (nx,ny,nz)
%       all_node: cell of nodes for each tissue layer's tetrahedral mesh, dimension (1,number of tissue layers)
%           each cell will have a node array containing the node coordinates of the mesh, dimension (nnodes,3)
%       all_elem: cell of elements for each tissue layer's tetrahedral mesh, dimension (1,number of tissue layers)
%           each cell will have an element array with the element list of the mesh, dimension (nnodes,5)
%       all_face: cell of faces for each tissue layer's tetrahedral mesh, dimension (1,number of tissue layers)
%           each cell will have an face array with the mesh surface element list of the tetrahedral mesh, dimension (nnodes,4)

%%

cd ..
addpath(genpath('.'))

volume_dir=['.' filesep 'volumes'];
toolbox_dir='.';
volume_name='stairhead.bin';
mesh_filename='stairhead_mesh.mat';

volume_structure=generate_mri_mesh(volume_dir,toolbox_dir,volume_name,mesh_filename);

%%