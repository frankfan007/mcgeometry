   function volume_structure=generate_mri_mesh(volume_direc,toolbox_direc,volume_name,mesh_filename,varargin);

%% generates mesh from volume file

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

%% checking if volume exists

if exist([volume_direc filesep volume_name],'file')
    direc=volume_direc;
    fileID=fopen([volume_direc filesep volume_name]);
elseif exist([toolbox_direc filesep volume_name],'file')
    direc=toolbox_direc;
    fileID=fopen([toolbox_direc filesep volume_name]);
else
    error('Volume not found in either volume directory or toolbox directory%s\n','');
end

%% reading volume

print_box('Reading volume...','',70);

read_vol=fread(fileID);
side_dim=nthroot(length(read_vol),3);

if side_dim==round(side_dim)
    permuted_vol=reshape(read_vol,[side_dim side_dim side_dim]);
else
    try
        permuted_vol=reshape(read_vol,[180 180 100]);
    catch
        prompt='Please input the dimensions of the volume in xyz coordinates: ';
        dlgtitle='Dimensions of volume';
        dims=[1 35];
        definput={'180 180 100'};
        opts.WindowStyle='normal';
        dim_vol=inputdlg(prompt,dlgtitle,dims,definput,opts);
        dim_vol=str2num(dim_vol{1});
        permuted_vol=reshape(read_vol,[dim_vol(1) dim_vol(2) dim_vol(3)]);
    end
end

% num_tissue_layers=length(unique(permuted_vol))-1;
%% creating and saving mesh

if exist([direc filesep mesh_filename],'file')
    subtext='loading mesh...';
    print_box('Mesh exists',subtext,70);
    load([direc filesep mesh_filename]);
else 
    
    subtext='creating mesh...';
    print_box('Mesh does not exist',subtext,70);
    
    tiss_indices=unique(permuted_vol);
    for tiss_type=transpose(tiss_indices(tiss_indices>0 & tiss_indices<6))
        [all_node{tiss_type},all_elem{tiss_type},all_face{tiss_type}]=v2m(permuted_vol,tiss_type,2,100);
    end
    
    node=all_node{1};
    elem=all_elem{1};
    face=all_face{1};
   
    save([direc filesep mesh_filename],'all_node','all_elem','all_face','node','elem','face','permuted_vol');
    
    if ~isempty(varargin)
        volume_name_vite=varargin{1};
        mesh_filename_vite=varargin{2};
        fileID=fopen([direc filesep volume_name_vite]);
        read_vol=fread(fileID);
        permuted_vol_vite=reshape(read_vol,[side_dim side_dim side_dim]);
        
        for tiss_type=transpose(tiss_indices(tiss_indices>0 & tiss_indices<6))
            [all_node_vite{tiss_type},all_elem_vite{tiss_type},all_face_vite{tiss_type}]=v2m(permuted_vol_vite,tiss_type,2,100);
        end
        
        node_vite=all_node_vite{1};
        elem_vite=all_elem_vite{1};
        face_vite=all_face_vite{1};
        
        volume_structure.vol_vite=permuted_vol_vite;
        
        volume_structure.all_node_vite=all_node_vite;
        volume_structure.all_elem_vite=all_elem_vite;
        volume_structure.all_face_vite=all_face_vite;
        
        save([direc filesep mesh_filename_vite],'all_node_vite','all_elem_vite','all_face_vite','node_vite','elem_vite','face_vite','permuted_vol_vite','-append');
    end
end

fclose all;

%% writing out variables

volume_structure.vol=permuted_vol;

volume_structure.all_node=all_node;
volume_structure.all_elem=all_elem;
volume_structure.all_face=all_face;

