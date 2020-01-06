function [node_closest_to_point,plane,point_unit_vec]=get_normal_plane_from_point(node,face,specified_point)
%
% [node_closest_to_point,plane,point_unit_vec]=get_normal_plane_from_point(node,face,specified_point)
%
% obtains normal plane to center of volume from specified point on volume surface

% input:
%   node: array containing node coordinates of mesh, dimension (nnodes,3)
%   face: array containing face coordinates of mesh, dimension (nnodes,4)
%   fiducial: array containing coordinates of specified point, dimension (1,3)

% output:
%   node_closest_to_point: coordinates of node closest to point, dimension (1,3)
%   plane: coordinates that determine normal plane from point, dimension (3,3)
%   point_unit_vec: unit normal vector from node closest to point

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

all_distances_from_isosurf=[];
for idx=1:size(node,1)
    all_distances_from_isosurf(idx)=norm(node(idx,:)-specified_point);
end

[~,sorted_indices]=sort(all_distances_from_isosurf);
node_closest_to_point=node(sorted_indices(1),:);
num_nodes_to_sample=20;

snorm=surfacenorm(node,face,'Normalize',1);

idx=1;
for single_node=1:num_nodes_to_sample
    I=sorted_indices(single_node);
    for face_column=1:3
        faces_with_source=find(face(:,face_column)==I);
        for f_index=1:length(faces_with_source)
            all_point_snorm(idx,:)=snorm(faces_with_source(f_index),:);
            idx=idx+1;
        end
    end
end

point_unit_vec=mean(all_point_snorm);

point1=node(sorted_indices(1),:);
point2=point1+point_unit_vec*10;
point3=[0 point1(2) point1(3)];

plane=[point1;point2;point3];
            
            