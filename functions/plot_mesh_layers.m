function plot_mesh_layers(all_node,all_face)
%
% plot_mesh_layers(all_node,all_face)
%
% plots mesh for volume with multiple tissue types
% input:
%   all_node: cell with arrays of nodes for each tissue layer, dimension (1, ntissue layers)
%   all_face: cell with arrays of faces for each tissue layer, dimension (1, ntissue layers)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

face_colors=[102 102 255;
    255 128 0;
    255 255 51;
    153 255 51;
    102 255 255;
    255 51 51;
    178 102 255;
    255 102 178;
    160 160 160]/256;

for tiss_type=1:length(all_node)
    trimesh(all_face{tiss_type}(:,1:3),all_node{tiss_type}(:,1),all_node{tiss_type}(:,2),all_node{tiss_type}(:,3),...
        'facecolor','none','EdgeColor',face_colors(tiss_type,:));
    hold on;
end
hold off
axis equal