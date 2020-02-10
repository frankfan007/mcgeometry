function generate_input_file(dir_struct,ref_param,mc_param)
% sets up input parameters and other variables from main script to call write_input_file function
% creates input file for Monte Carlo simulation
%
% generate_input_file(dir_struct,ref_param,mc_param)
%
% input:
%   dir_struct: structure with fields
%       volume_dir: directory where volume is stored
%       dcs_mc_toolbox: directory where DCS and MC functions and scripts are stored
%       input_filename: full filename of input file
%   mc_param: structure with fields
%       mus: array with scattering values for each tissue layer, dimension (1, number of tissue layers)
%       volume_name: name of volume,
%   ref_param: substructure with fields
%       all_locations: array with the source position as the first row, and the detector positions as the subsequent rows
%           dimension (1+number of detectors, 3)
%       source_unit_vec: array with the unit normal vector from the source position, dimension (1,3)
%       fiducial_point: fiducial point position, dimension (1,3)
%       vol: volume array, dimension (nx,ny,nz)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%%
% setting variables
volume_direc=dir_struct.volume_dir;
toolbox_direc=[dir_struct.dcs_mc_toolbox filesep 'volumes'];

% extracting number of tissue types 
if exist([volume_direc filesep mc_param.volume_name],'file')
    direc=volume_direc;
    mc_param.volume_name_unix=[dir_struct.volume_dir_unix filesep mc_param.volume_name];
    mc_param.volume_name_unix=strrep(mc_param.volume_name_unix,'\','/');
    fileID=fopen([volume_direc filesep mc_param.volume_name]);
elseif exist([toolbox_direc filesep mc_param.volume_name],'file')
    direc=toolbox_direc;
    mc_param.volume_name_unix=[dir_struct.dcs_mc_toolbox_unix 'volumes/' mc_param.volume_name];
    fileID=fopen([toolbox_direc filesep mc_param.volume_name]);
end

file_contents=fread(fileID);
mc_param.num_tissue_layers=max(unique(file_contents));
fclose(fileID);

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

