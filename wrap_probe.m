function ref_param=wrap_probe(ref_param)
%
% ref_param=wrap_probe(ref_param)
%
% determines source and detector locations on a given volume
%
% input:
%   ref_param: structure with fields
%       vol: volume
%       all_node: 1 x n cell of nodes for each tissue layer, where n is the number of tissue layers
%       all_elem: 1 x n cell of elements for each tissue layer, where n is the number of tissue layers
%       all_face: 1 x n cell of faces for each tissue layer, where n is the number of tissue layers
%       det_distances: 1 x d array of distances in mm for detectors, where d is number of detectors
%       default_fiducial_pos: 1 x 3 source position in generic onionhead
%
% output:
%   ref_param: structure with added fields
%       outer_line: 1 x r cell of outer loops that the source/detectors will lie on,
%           where r is the number of rotations
%       source: r x 3 array of source locations, where r is the number of rotations
%       source_unit_vec: r x 3 array of source unit vectors, where r is the number of rotations
%       det_arr: r x d x 3 array of detector locations, where r is the number of rotations and
%           d is the number of detectors
%       det_unit_vec: r x d x 3 array of detector unit vectors, where r is the number of rotations and
%           d is the number of detectors
%       all_locations: r x (d+1) x 3 array of source and detector locations, where the first position of
%           the second index is the source
%       fiducial_point: fiducial location (manually chosen or fiducial driven)
%       rotate_deg: array of rotation degrees

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%% setting variables

fv=isosurface(ref_param.vol,0.5,'noshare');

node=ref_param.all_node;
elem=ref_param.all_elem;
face=ref_param.all_face;

scalp_node=node{1};
scalp_elem=elem{1};
scalp_face=face{1};

if isfield(ref_param,'all_node_vite')
    node_vite=ref_param.all_node_vite;
    elem_vite=ref_param.all_elem_vite;
    face_vite=ref_param.all_face_vite;
    
    scalp_node_vite=node_vite{1};
    scalp_elem_vite=elem_vite{1};
    scalp_face_vite=face_vite{1};
end

snorm=surfacenorm(scalp_node,scalp_face,'Normalize',1);
%% determine fiducial location

dlgtitle = 'Input';
definput={'180 80 200'};
opts.WindowStyle='normal';
dims = [1 35];

if ref_param.use_default_fiducial
    fiducial=ref_param.default_fiducial_pos;
    [node_closest_to_fiducial,plane,fiducial_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,fiducial);
    print_box('Using default fiducial position','',70);
else
    if isempty(ref_param.fiducial_pos)
        if isfield(ref_param,'all_node_vite')
            idx=1;
            while true
                vite_mesh=figure(305);
                trimesh(scalp_face_vite(:,1:3),scalp_node_vite(:,1),scalp_node_vite(:,2),scalp_node_vite(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
                axis equal
                rotate3d on
                hold on
                if idx==1 && exist('node_closest_to_default','var')
                    plot3(node_closest_to_default(1),node_closest_to_default(2),node_closest_to_default(3),'ro','LineWidth',2,'MarkerFaceColor','r')
                end
                hold off
                
                dcm_obj=datacursormode(vite_mesh);
                set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');
        
                vite_mesh.CurrentCharacter='p';
        
                m=msgbox('Please select data point on the mesh with data point cursor, then press space when satisfied');
        
                waitfor(vite_mesh,'CurrentCharacter',char(32));
                dcm_obj=datacursormode(vite_mesh);
                cursor_info=getCursorInfo(dcm_obj);
                
                if length(cursor_info)>1
                    too_many_points=msgbox('Please only select one fiducial location');
                    close(too_many_points);
                    continue
                elseif isempty(cursor_info)
                    no_points_selected=msgbox('Please select at least one fiducial location');
                    close(no_points_selected);
                    continue
                end
                
                fiducial=cursor_info.Position;
                
                [node_closest_to_fiducial,plane,fiducial_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,fiducial);
                [cutpos_vite,cutvalue_vite,facedata_vite]=qmeshcut(scalp_elem_vite(:,1:4),scalp_node_vite,scalp_node_vite(:,1),plane);
                
                figure(305);
                trimesh(scalp_face_vite(:,1:3),scalp_node_vite(:,1),scalp_node_vite(:,2),scalp_node_vite(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
                hold on;
                plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'ro','LineWidth',2,'MarkerFaceColor','r')
                patch('Faces',facedata_vite,'Vertices',cutpos_vite,'FaceVertexCData',cutvalue_vite,'facecolor','k');
                axis equal
                view(-180,90); rotate3d on
                hold off
                
                satisfied_dlg=inputdlg('Are you satisfied with this cut?','Satisfied with cut?',dims,{'yes'},opts);
                if isempty(satisfied_dlg)
                    return
                elseif strmp(satisfied_dlg{1},'yes')
                    break
                end
                idx=idx+1;
            end
        else
            surfnodes=unique(reshape(scalp_face(:,1:3),[size(scalp_face,1)*3 1]));
            scalp_surf=scalp_node(surfnodes,:);
            
            all_distances_from_default=[];
            for idx=1:size(scalp_surf,1)
                all_distances_from_default(idx)=norm(scalp_surf(idx,:)-ref_param.default_fiducial_pos);
            end
            
            [~,I]=min(all_distances_from_default);
            node_closest_to_default=scalp_surf(I,:);
            
            figure(300);
            trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
            axis equal
            rotate3d on
            hold on;
            plot3(node_closest_to_default(1),node_closest_to_default(2),node_closest_to_default(3),'ro','LineWidth',2,'MarkerFaceColor','r')
            
            use_default=inputdlg('Use default fiducial location? Type "yes" for yes and "no" for no', 'Default location',dims,{'yes'},opts);
            
            if isempty(use_default)
                return
            elseif strcmp(use_default,'yes')
                fiducial=node_closest_to_default;
                [node_closest_to_fiducial,plane,fiducial_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,fiducial);
                [cutpos,cutvalue,facedata]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),plane);
                
                figure(300)
                patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k');
            else
                idx=1;
                while true
                    meshfig=figure(300);
                    trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
                    axis equal
                    rotate3d on
                    hold on;
                    if idx==1 && exist('node_closest_to_default','var')
                        plot3(node_closest_to_default(1),node_closest_to_default(2),node_closest_to_default(3),'ro','LineWidth',2,'MarkerFaceColor','r')
                    end
                    hold off
                    
                    dcm_obj=datacursormode(meshfig);
                    set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on');
        
                    meshfig.CurrentCharacter='p';
        
                    m=msgbox('Please select data point on the mesh with data point cursor, then press space when satisfied');
        
                    waitfor(meshfig,'CurrentCharacter',char(32));
                    dcm_obj=datacursormode(meshfig);
                    cursor_info=getCursorInfo(dcm_obj);
                    fiducial=cursor_info.Position;
                    
                    [node_closest_to_fiducial,plane,fiducial_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,fiducial);
                    [cutpos,cutvalue,facedata]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),plane);
                    
                    figure(300);
                    trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
                    axis equal; rotate3d on; hold on; view(-180,90);
                    plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'bo','LineWidth',2,'MarkerFaceColor','b')
                    patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k');
                    hold off
                    
                    satisfied_dlg=inputdlg('Are you satisfied with this cut?','Satisfied with cut?',dims,{'yes'},opts);
                    if isempty(satisfied_dlg)
                        return
                    elseif strcmp(satisfied_dlg{1},'yes')
                        break
                    end
                    idx=idx+1;
                end
            end
        end
    else
        fiducial=ref_param.fiducial_pos;
        [node_closest_to_fiducial,plane,fiducial_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,fiducial);
        print_box('Using specified fiducial position','',70);
    end
end

[cutpos,cutvalue,facedata]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),plane);

 %% makes plane from fiducial point chosen above
 
 if isempty(ref_param.rotate_deg)
     while true
         figure(100);
         trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
         hold on;
         plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'ro','LineWidth',2,'MarkerFaceColor','r')
         hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k');
         axis equal;
         view(-180,90); rotate3d on
         
         if isfield(ref_param,'all_node_vite')
             [cutpos_vite,cutvalue_vite,facedata_vite]=qmeshcut(scalp_elem_vite(:,1:4),scalp_node_vite,scalp_node_vite(:,1),plane);
             
             figure(105);
             trimesh(scalp_face_vite(:,1:3),scalp_node_vite(:,1),scalp_node_vite(:,2),scalp_node_vite(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
             hold on;
             plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'ro','LineWidth',2,'MarkerFaceColor','r')
             hcut_vite=patch('Faces',facedata_vite,'Vertices',cutpos_vite,'FaceVertexCData',cutvalue_vite,'facecolor','k');
             axis equal
             view(-180,90); rotate3d on
             hold off
         end
         
         definput={'0'};
         rot_ans=inputdlg('Enter rotation degree: ',dlgtitle,dims,definput,opts);
         deg=str2num(rot_ans{1});
         
         figure(100);
         hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k','visible','off');
         rotate(hcut,fiducial_unit_vec,deg,node_closest_to_fiducial);
         
         newplane=hcut.Vertices(1:3,:);
         [newcutpos,newcutvalue,newfacedata]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),newplane);
         hcut_new=patch('Faces',newfacedata,'Vertices',newcutpos,'FaceVertexCData',newcutvalue,'facecolor','k');
         hold off
         
         if isfield(ref_param,'all_node_vite')
             
             figure(105);
             hcut_vite=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k','visible','off');
             rotate(hcut_vite,fiducial_unit_vec,deg,node_closest_to_fiducial);
             
             newplane_vite=hcut_vite.Vertices(1:3,:);
             [newcutpos_vite,newcutvalue_vite,newfacedata_vite]=qmeshcut(scalp_elem_vite(:,1:4),scalp_node_vite,scalp_node_vite(:,1),newplane_vite);
             hcut_new_vite=patch('Faces',newfacedata_vite,'Vertices',newcutpos_vite,'FaceVertexCData',newcutvalue_vite,'facecolor','k');
             hold off
             
         end
         
         satisfied_dlg=inputdlg('Are you satisfied with this rotation?','Satisfied with cut?',dims,{'yes'},opts);
         if isempty(satisfied_dlg)
            return
         elseif strcmp(satisfied_dlg{1},'yes')
             break
         end
     end
 else
     deg=ref_param.rotate_deg;
     
     figure(100);
     trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
     hold on;
     plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'ro','LineWidth',2,'MarkerFaceColor','r')
     hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k');
     axis equal;
     view(-180,90); rotate3d on
     
     hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k','visible','off');
     rotate(hcut,fiducial_unit_vec,deg,node_closest_to_fiducial);
     
     newplane=hcut.Vertices(1:3,:);
     [newcutpos,newcutvalue,newfacedata]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),newplane);
     hcut_new=patch('Faces',newfacedata,'Vertices',newcutpos,'FaceVertexCData',newcutvalue,'facecolor','k');
     hold off
     
     if isfield(ref_param,'all_node_vite')
         [cutpos_vite,cutvalue_vite,facedata_vite]=qmeshcut(scalp_elem_vite(:,1:4),scalp_node_vite,scalp_node_vite(:,1),plane);
         
         figure(105);
         trimesh(scalp_face_vite(:,1:3),scalp_node_vite(:,1),scalp_node_vite(:,2),scalp_node_vite(:,3),'facecolor','none','EdgeColor',[255 128 0]/256)
         hold on;
         plot3(node_closest_to_fiducial(1),node_closest_to_fiducial(2),node_closest_to_fiducial(3),'ro','LineWidth',2,'MarkerFaceColor','r')
         hcut_vite=patch('Faces',facedata_vite,'Vertices',cutpos_vite,'FaceVertexCData',cutvalue_vite,'facecolor','k');
         axis equal
         view(-180,90); rotate3d on
         hold off
         
         hcut_vite=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','k','visible','off');
         rotate(hcut_vite,fiducial_unit_vec,deg,node_closest_to_fiducial);
         
         newplane_vite=hcut_vite.Vertices(1:3,:);
         [newcutpos_vite,newcutvalue_vite,newfacedata_vite]=qmeshcut(scalp_elem(:,1:4),scalp_node,scalp_node(:,1),newplane_vite);
         hcut_new_vite=patch('Faces',newfacedata_vite,'Vertices',newcutpos_vite,'FaceVertexCData',newcutvalue_vite,'facecolor','k');
         hold off
     end
 end


%% make outline of plane

print_box('Placing probe on head...','finding orientation and source-detector locations',70);

[bcutpos,~,bcutedges]=qmeshcut(scalp_face(:,1:3),scalp_node,scalp_node(:,1),newplane);
[bcutpos,bcutedges]=removedupnodes(bcutpos,bcutedges);
bcutloop=extractloops(bcutedges);

bcutloop(isnan(bcutloop))=[]; % there can be multiple loops, remove the separators

outer_line=bcutpos(bcutloop,:);

for point=1:size(outer_line,1)
    all_distances_outer_line(point)=norm(fiducial-outer_line(point,:));
end

[~,index]=min(all_distances_outer_line);

%% adjust starting point of the outer line until the fiducial is not close to the start

shift_tolerance=50;
shift_length=100;
interv=20;

% shift if too close to fiducial
if index<shift_tolerance || (size(outer_line,1)-index)<shift_tolerance
    outer_line_shift=outer_line(shift_length:end,:);
    outer_line_shift=cat(1,outer_line_shift,outer_line(1:(shift_length-1),:));
    for point=1:size(outer_line_shift,1)
        all_distances_outer_line_shift(point)=norm(fiducial-outer_line_shift(point,:));
    end
    [~,index_shift]=min(all_distances_outer_line_shift);
else
    outer_line_shift=outer_line;
    index_shift=index;
end

count=0;
err_count=0;
while count==err_count
    try
        [x,xa,~]=unique(outer_line_shift((index_shift-interv):(index_shift+interv),1));
        y=outer_line_shift((index_shift-interv):(index_shift+interv),2);
        y=y(xa);
        z=outer_line_shift((index_shift-interv):(index_shift+interv),3);
        z=z(xa);
    catch
        shift_length=shift_length+10;
        interv=interv+5;
        
        outer_line_shift=outer_line(shift_length:end,:);
        outer_line_shift=cat(1,outer_line_shift,outer_line(1:(shift_length-1),:));
        for point=1:size(outer_line_shift,1)
            all_distances_outer_line_shift(point)=norm(fiducial-outer_line_shift(point,:));
        end
            [~,index_shift]=min(all_distances_outer_line_shift);
        
        err_count=err_count+1;
    end
    count=count+1;
end

%%  figure out orientation of head

minx=min(scalp_node(:,1)); maxx=max(scalp_node(:,1));
miny=min(scalp_node(:,2)); maxy=max(scalp_node(:,2));
minz=min(scalp_node(:,3)); maxz=max(scalp_node(:,3));

if (maxx-minx)<(maxy-miny)
    idx_to_compare=1;
    midpoint=minx+((maxx-minx)/2);
else
    idx_to_compare=2;
    midpoint=miny+((maxy-miny)/2);
end

%% take section of the outer line surrounding the fiducial, and interpolate then smooth it

orig_points=[x,y,z];

CS = cat(1,0,cumsum(sqrt(sum(diff(orig_points,[],1).^2,2))));
interp_points_jagged = interp1(CS, orig_points, unique([CS(:)' linspace(0,CS(end),1000)]),'spline');

[interp_points,~]=smoothdata(interp_points_jagged,1,'sgolay');

arr_to_cmp=interp_points(:,idx_to_compare);

all_differences=midpoint-arr_to_cmp;

if (sum(all_differences>0)/length(all_differences)>0.5 && ismonotonic(arr_to_cmp,[],'decreasing')) || (sum(all_differences<0)/length(all_differences)>0.5 && ismonotonic(arr_to_cmp,[],'increasing'))
    interp_points=flipud(interp_points);
end

outline=interp_points;

%% find source and detector locations

sd_sep=ref_param.det_distances(end)/2;

for pt=1:size(interp_points,1)
    distance_from_fiducial(pt)=norm(node_closest_to_fiducial-interp_points(pt,:));
end

source_start_idx=find(abs(distance_from_fiducial-sd_sep)<0.1,1);
source=interp_points(source_start_idx,:);
temp_interp_points=interp_points(source_start_idx:end,:);

for det=2:size(temp_interp_points,1)
    distance_between_pts(det-1)=norm(temp_interp_points(det,:) - temp_interp_points(det-1,:));
end

sum_dist=cumsum(distance_between_pts);

for idx=1:length(ref_param.det_distances)
    [~,index_values(idx)]=min(abs(sum_dist-ref_param.det_distances(idx)));
end

det_arr=temp_interp_points(index_values,:);

%% combining into single array and finding normals for source and each detector

locations(1,:)=source;

% find source normal
[~,~,source_unit_vec]=get_normal_plane_from_point(scalp_node,scalp_face,source);

% find detector normals
for det_idx=1:length(ref_param.det_distances)
    locations(det_idx+1,:)=det_arr(det_idx,:);
    [~,~,det_unit_vec(det_idx,:)]=get_normal_plane_from_point(scalp_node,scalp_face,squeeze(det_arr(det_idx,:)));
end

%% plotting result

figure(200)
trimesh(scalp_face(:,1:3),scalp_node(:,1),scalp_node(:,2),scalp_node(:,3),'facecolor','none')
hold on; axis equal;
plot3(fiducial(1),fiducial(2),fiducial(3),'ko','LineWidth',2,'MarkerFaceColor','k');
plot3(interp_points(:,1),interp_points(:,2),interp_points(:,3),'k-','LineWidth',1);
plot3(source(1),source(2),source(3),'ro','LineWidth',2,'MarkerFaceColor','r');
plot3(squeeze(det_arr(:,1)),squeeze(det_arr(:,2)),squeeze(det_arr(:,3)),'bo','LineWidth',2,'MarkerFaceColor','b');
quiver3(source(1),source(2),source(3),source_unit_vec(1)*50,source_unit_vec(2)*50,source_unit_vec(3)*50,'LineWidth',2,'Color','r');
quiver3(det_arr(:,1),det_arr(:,2),det_arr(:,3),...
    squeeze(det_unit_vec(:,1))*70,squeeze(det_unit_vec(:,2))*70,squeeze(det_unit_vec(:,3))*70,'LineWidth',2,'Color','b');
drawnow

if isfield(ref_param,'all_node_vite')
    figure(205);
    trimesh(scalp_face_vite(:,1:3),scalp_node_vite(:,1),scalp_node_vite(:,2),scalp_node_vite(:,3),'facecolor','none')
    hold on
    axis equal
    plot3(fiducial(1),fiducial(2),fiducial(3),'ko','LineWidth',2,'MarkerFaceColor','k');
    view(-180,90); rotate3d on
    hold on;
    plot3(interp_points(:,1),interp_points(:,2),interp_points(:,3),'k-','LineWidth',1);
    plot3(source(1),source(2),source(3),'ro','LineWidth',2,'MarkerFaceColor','r');
    plot3(squeeze(det_arr(:,1)),squeeze(det_arr(:,2)),squeeze(det_arr(:,3)),'bo','LineWidth',2,'MarkerFaceColor','b');
    quiver3(source(1),source(2),source(3),source_unit_vec(1)*50,source_unit_vec(2)*50,source_unit_vec(3)*50,'LineWidth',2,'Color','r');
    quiver3(det_arr(:,1),det_arr(:,2),det_arr(:,3),...
        squeeze(det_unit_vec(:,1))*70,squeeze(det_unit_vec(:,2))*70,squeeze(det_unit_vec(:,3))*70,'LineWidth',2,'Color','b');
    axis equal
    view(-180,90); rotate3d on
    hold off
    drawnow
end

%% moving points so that they are on isosurface

for pt=1:size(locations,1)
    single_point=squeeze(locations(pt,:));
    single_point=[single_point(2) single_point(1) single_point(3)];
    for vert_idx=1:size(fv.vertices,1)
        all_distances_from_isosurf(vert_idx)=norm(fv.vertices(vert_idx,:)-single_point);
    end
    [all_mins(pt),all_indices(pt)]=min(all_distances_from_isosurf);
    clear all_distances_from_isosurf
end

for pt=1:size(all_indices,2)
    adj_locations(pt,:)=[fv.vertices(all_indices(pt),2) fv.vertices(all_indices(pt),1) fv.vertices(all_indices(pt),3)];
end

%% plotting isosurface

figure
isosurface(ref_param.vol,0.5);
axis equal; grid on; hold on;
plot3(squeeze(adj_locations(:,2)),squeeze(adj_locations(:,1)),squeeze(adj_locations(:,3)),'k*');
drawnow

%% saving variables in a reference parameter structure

ref_param.outer_line=outline;
ref_param.source=source;
ref_param.source_unit_vec=source_unit_vec;
ref_param.det_unit_vec=det_unit_vec;
ref_param.det_arr=det_arr;
ref_param.all_locations=adj_locations;
ref_param.fiducial_point=fiducial;
ref_param.rotate_deg=deg;

%%
