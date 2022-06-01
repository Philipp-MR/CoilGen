function coil_parts= create_sweep_along_surface(coil_parts,input)
%create a volumetric coil body by surface sweep, @Philipp Amrein, Uniklinik
%Freiburg 2022


for part_ind=1:numel(coil_parts)

parameterized_mesh=coil_parts(part_ind).coil_mesh;
points_to_shift=coil_parts(part_ind).points_to_shift;
shift_array=coil_parts(part_ind).shift_array;
wire_path=coil_parts(part_ind).wire_path;
planary_mesh_matlab_format=triangulation(parameterized_mesh.faces',parameterized_mesh.uv');
curved_mesh_matlab_format=triangulation(parameterized_mesh.faces',parameterized_mesh.v);
cross_section_points=input.cross_sectional_points;
save_mesh=input.save_stl_flag;
output_directory=input.output_directory;
conductor_conductivity=input.specific_conductivity_conductor;

convolutional_vector_length=1; %for smothering the curverture along the track

%center the cross section around the [0,0] origin that roation will be
%later valid
cross_section_points=cross_section_points-mean(cross_section_points,2);
%build a triangulation from the cross section
cross_section_points=cross_section_points(:,1:end-1);

%build a 2d mesh of the cross section by the corner points
num_cross_section_points=size(cross_section_points,2);
cross_section_edges=[1:num_cross_section_points; circshift(1:num_cross_section_points,-1)];
cross_section_triangulation=delaunayTriangulation(cross_section_points',cross_section_edges');
cross_section_points=[cross_section_triangulation.Points'; zeros(1,size(cross_section_triangulation.Points,1))]; %embed the cross section in 3d
is_interior_tri = isInterior(cross_section_triangulation);
cross_section_triangulation=triangulation(cross_section_triangulation.ConnectivityList(is_interior_tri,:),cross_section_triangulation.Points);

%Calaculate the area of the cross section
cross_section_triangulation_3d=triangulation(cross_section_triangulation.ConnectivityList,[cross_section_triangulation.Points zeros(size(cross_section_triangulation.Points,1),1)]);
cross_section_area=0;
for tri_ind=1:size(cross_section_triangulation.ConnectivityList,1)
    P1=cross_section_triangulation_3d.Points(cross_section_triangulation_3d.ConnectivityList(tri_ind,1),:)';
    P2=cross_section_triangulation_3d.Points(cross_section_triangulation_3d.ConnectivityList(tri_ind,2),:)';
    P3=cross_section_triangulation_3d.Points(cross_section_triangulation_3d.ConnectivityList(tri_ind,3),:)';
    cross_section_area = cross_section_area+1/2*norm(cross(P2-P1,P3-P1));
end
%Calculate the length of the coil
wire_path.v_length=sum(vecnorm(wire_path.v(:,2:end)-wire_path.v(:,1:end-1)));
%Calculate the ohmian resistance
ohmian_resistance=wire_path.v_length/(cross_section_area*conductor_conductivity);

%calculate a radius of the conductor cross section which is later
%important to avoid intersection between angulated faces
cross_section_center=sum(cross_section_points,2)./size(cross_section_points,2);
cross_section_radius=max(vecnorm(cross_section_points-cross_section_center));

%Remove repeating entries
wire_path.v(:,find(vecnorm(diff(wire_path.v,1,2))==0))=[]; 
wire_path.uv(:,find(vecnorm(diff(wire_path.uv,1,2))==0))=[]; 

%Open the track if its not already opened
point_inds_to_delete=vecnorm(wire_path.v(:,end)-wire_path.v)<cross_section_radius/2;
point_inds_to_delete=find(point_inds_to_delete(1:round(numel(point_inds_to_delete)/2)));
wire_path.v(:,find(point_inds_to_delete))=[];
wire_path.v(:,find(point_inds_to_delete))=[];

% if all(wire_path.v(:,1)==wire_path.v(:,end))
%     wire_path.v(:,end)=[];
%     wire_path.uv(:,end)=[];
% end
% %Close the track if its not already closed
% if all(wire_path.v(:,1)~=wire_path.v(:,end))
%     wire_path.v=[wire_path.v wire_path.v(:,1)];
%     wire_path.uv=[wire_path.uv wire_path.uv(:,1)];
% end



% cross_section_points_rep=repmat(cross_section_points,[1,size(wire_path.v,2)]);
% wire_path.v_rep=repelem(wire_path.v,1,size(cross_section_points,2));
% sweeped_face=cross_section_points_rep+wire_path.v_rep;


% figure; hold on; axis equal;
% for sweep_ind=1:size(wire_path.v,2)-1
% plot_inds=((sweep_ind-1)*num_cr_points+1):(sweep_ind)*num_cr_points;
% plot3(sweeped_face(1,plot_inds),sweeped_face(2,plot_inds),sweeped_face(3,plot_inds),'r');
% end


%calculate the normal vectors along the wire track
surface_normal_alonge_wire_path.v=zeros(3,size(wire_path.v,2));
for point_ind=1:size(wire_path.v,2)
node_ind_normals_target = pointLocation(planary_mesh_matlab_format,wire_path.uv(1,point_ind),wire_path.uv(2,point_ind));
if isnan(node_ind_normals_target) %make excepetions for strange output of pointLocation
surface_normal_alonge_wire_path.v(:,point_ind)=surface_normal_alonge_wire_path.v(:,point_ind-1);
else
surface_normal_alonge_wire_path.v(:,point_ind)=parameterized_mesh.fn(node_ind_normals_target,:)';
end
end
% smooth the normals that there are no sharp twists
conv_vec=[0:convolutional_vector_length convolutional_vector_length-1:-1:0]./convolutional_vector_length;
% surface_normal_alonge_wire_path.v = conv2(surface_normal_alonge_wire_path.v,conv_vec,'same');
% surface_normal_alonge_wire_path.v=surface_normal_alonge_wire_path.v./repmat(vecnorm(surface_normal_alonge_wire_path.v),[3 1]); % normalize 




% %calculate the normal vectors
% surface_normal_alonge_wire_path.v=zeros(3,size(wire_path.v,2));
% for point_ind=1:size(wire_path.v,2)
% node_ind_normals_target = nearestNeighbor(curved_mesh_matlab_format,wire_path.v(1,point_ind),wire_path.v(2,point_ind),wire_path.v(3,point_ind)); 
% surface_normal_alonge_wire_path.v(:,point_ind)=parameterized_mesh.n(node_ind_normals_target,:)';
% end
% % smooth the normals that there are no sharp twists
% conv_vec=[0:convolutional_vector_length convolutional_vector_length-1:-1:0]./convolutional_vector_length;
% surface_normal_alonge_wire_path.v = conv2(surface_normal_alonge_wire_path.v,conv_vec,'same');
% surface_normal_alonge_wire_path.v=surface_normal_alonge_wire_path.v./repmat(vecnorm(surface_normal_alonge_wire_path.v),[3 1]); % normalize 



%prepare the track direction
path_directions=wire_path.v(:,2:end)-wire_path.v(:,1:end-1); %edge vectors
path_directions=[path_directions path_directions(:,end)]; %add a repetition at the end for the last face
path_directions=path_directions./repmat(vecnorm(path_directions),[3 1]); % normalize

path_directions(:,2:end-1)=(path_directions(:,1:end-2)+path_directions(:,3:end))./2; % average the vector between its precessor and suczsecor
path_directions=path_directions./repmat(vecnorm(path_directions),[3 1]); % normalize again

%smooth the path that there are no sharp turns
% path_directions = conv2(path_directions,conv_vec,'same');
% path_directions=path_directions./repmat(vecnorm(path_directions),[3 1]); % normalize 

%%%%%sweeping of the surface along the path
%defining the face normal of the start of the conductor
face_points=zeros([size(wire_path.v,2),size(wire_path.v,1),size(cross_section_points,2)]); %initialize the face points
all_node_points=zeros(size(wire_path.v,1)*size(cross_section_points,2),3);
run_ind=1;
for point_ind=1:size(wire_path.v,2)
        e1=path_directions(:,point_ind);
        e2=surface_normal_alonge_wire_path.v(:,point_ind);
        e3=cross(e1,e2)./norm(cross(e1,e2));
        e2=cross(e1,e3)./norm(cross(e1,e3));
        %Rotate the face points to have the right orientation
        for aaaa=1:size(face_points,3)
        face_points(point_ind,:,aaaa)=[e2 e3 e1]*cross_section_points(:,aaaa);
        end
        %shift the oriented corner points to the position in the wire path
        for aaaa=1:size(face_points,3)
        face_points(point_ind,:,aaaa)=face_points(point_ind,:,aaaa)+wire_path.v(:,point_ind)';
        all_node_points(run_ind,:)=face_points(point_ind,:,aaaa);
        run_ind=run_ind+1;
        end
        
end


%Build the 3D mesh by sweeping the surface by only using the outer edges

%buidling the shell triangles that form the surface of the sweeped body
sweeped_surface_triangles=zeros((size(face_points,1)-1)*size(face_points,3)*2,3);
sweeped_surface_vertices=all_node_points;
full_edge_inds=[1:size(face_points,3) 1];
full_track_inds=[1:size(face_points,1)-1 1];
num_corners=size(face_points,3);

run_ind=2;
for track_ind=0:(size(face_points,1)-2)
    for edge_ind=1:size(face_points,3)
    node_a=track_ind*num_corners+full_edge_inds(edge_ind);
    node_b=track_ind*num_corners+full_edge_inds(edge_ind+1);
    node_c=full_track_inds((track_ind+1))*num_corners+full_edge_inds(edge_ind);
    node_d=full_track_inds((track_ind+1))*num_corners+full_edge_inds(edge_ind+1);
    tri_1=[node_a node_b node_d];
    tri_2=[node_d node_c node_a];   
    sweeped_surface_triangles(run_ind-1,:)=tri_1;
    sweeped_surface_triangles(run_ind,:)=tri_2;
    run_ind=run_ind+2;
    end
end


%Build the final triangles and close the surface
run_ind=2;
final_triangles=zeros(2*num_corners,3);
left_end_nodes=((size(wire_path.v,2)-1)*num_corners+1):(size(wire_path.v,2)*num_corners);
right_end_nodes=1:num_corners;

left_end_nodes=[left_end_nodes left_end_nodes(1)];
right_end_nodes=[right_end_nodes right_end_nodes(1)];

for edge_ind=1:num_corners
    node_a=left_end_nodes(edge_ind);
    node_b=left_end_nodes(edge_ind+1);
    node_c=right_end_nodes(edge_ind);
    node_d=right_end_nodes(edge_ind+1);
    tri_1=[node_a node_b node_d];
    tri_2=[node_d node_c node_a];
    final_triangles(run_ind-1,:)=tri_1;
    final_triangles(run_ind,:)=tri_2;
    run_ind=run_ind+2;
end

 layout_surface_mesh=triangulation([sweeped_surface_triangles; final_triangles],sweeped_surface_vertices);        

% % % % %building the mesh from the corner points
% % % % num_tri_per_face=size(cross_section_triangulation.ConnectivityList,1);
% % % % num_cross_face_nodes=size(face_points,3);
% % % % sweeped_tetras=zeros((size(face_points,1)-1)*num_tri_per_face*3,4);
% % % % sweeped_corners=zeros(3,size(face_points,1)*num_cross_face_nodes);
% % % % 
% % % % for face_ind=1:size(face_points,1)-1
% % % % for tri_ind=1:num_tri_per_face
% % % % 
% % % % %for each triangle in the cross section build a triangular prism to the
% % % % %next face which is composed out of 3 tetrahedral elements
% % % % 
% % % % cross_tri_inds=cross_section_triangulation.ConnectivityList +(face_ind-1)*num_cross_face_nodes;
% % % % 
% % % % %triangle corner inds of first triangle
% % % % ind_a=cross_tri_inds(tri_ind,1);
% % % % ind_b=cross_tri_inds(tri_ind,2);
% % % % ind_c= cross_tri_inds(tri_ind,3);
% % % % 
% % % % %triangle corner inds of second triangle
% % % % ind_d=cross_tri_inds(tri_ind,1)+num_cross_face_nodes;
% % % % ind_e=cross_tri_inds(tri_ind,2)+num_cross_face_nodes;
% % % % ind_f=cross_tri_inds(tri_ind,3)+num_cross_face_nodes;
% % % % 
% % % % sweeped_tetras((face_ind-1)*num_tri_per_face*3+(tri_ind-1)*3+1,:)=[ind_a ind_b ind_c ind_d];
% % % % sweeped_tetras((face_ind-1)*num_tri_per_face*3+(tri_ind-1)*3+2,:)=[ind_b ind_c ind_d ind_e];
% % % % sweeped_tetras((face_ind-1)*num_tri_per_face*3+(tri_ind-1)*3+3,:)=[ind_c ind_d ind_e ind_f];
% % % % 
% % % % end
% % % % end
% % % % 
% % % % 
% % % % %assign the corner coords to the sweeped body
% % % % for face_ind=1:size(face_points,1)
% % % %     sweeped_corners(:,((face_ind-1)*num_cross_face_nodes+1):((face_ind)*num_cross_face_nodes))=squeeze(face_points(face_ind,:,:));
% % % % end
% % % % 
% % % % %layout_volume_mesh=triangulation(sweeped_tetras,sweeped_corners');
% % % % 
% % % % %convert the volume tetrahedral mesh into a surface mesh
% % % % %surface_triangles = [sweeped_tetras(:,[1 2 3]);sweeped_tetras(:,[1 2 4]);sweeped_tetras(:,[1 3 4]);sweeped_tetras(:,[2 3 4])];
% % % % surface_triangles = [sweeped_tetras(:,[2 3 4]);sweeped_tetras(:,[1 4 3]);sweeped_tetras(:,[1 2 4]);sweeped_tetras(:,[1 3 2])];
% % % % 
% % % % %delete the inner faces (non-unique triangles)
% % % % sorted_surface_triangles = sort(surface_triangles,2);
% % % % [sorted_surface_triangles,sort_inds] = sortrows(sorted_surface_triangles);
% % % % surface_triangles=surface_triangles(sort_inds,:);
% % % % nun_unique_tri = find(all(diff(sorted_surface_triangles) == 0,2));
% % % % surface_triangles([nun_unique_tri;nun_unique_tri + 1],:) = [];
% % % % surface_triangles=surface_triangles(:,[1 3 2]);
% % % % 
% % % % layout_surface_mesh=triangulation(surface_triangles,sweeped_corners');        


%Save the mesh as an .stl file
if save_mesh
stlwrite(layout_surface_mesh,strcat(output_directory,'\sweeped_layout_part',num2str(part_ind),'.stl'),'text') 
stlwrite(curved_mesh_matlab_format,strcat(output_directory,'\surface_part',num2str(part_ind),'.stl'),'text') 
end

%Assign outputs
coil_parts(part_ind).layout_surface_mesh=layout_surface_mesh;
coil_parts(part_ind).ohmian_resistance=ohmian_resistance;


% figure;
% hold on;
% axis equal;
% %trisurf(triangulation(parameterized_mesh.f,parameterized_mesh.v),'facecolor','cyan','facealpha',0.5);
% trisurf(layout_surface_mesh,'facecolor','red') ; 
% % plot3(wire_path.v(1,:),wire_path.v(2,:),wire_path.v(3,:),'r','LineWidth',2);
% %quiver3(wire_path.v(1,:),wire_path.v(2,:),wire_path.v(3,:),surface_normal_alonge_wire_path.v(1,:),surface_normal_alonge_wire_path.v(2,:),surface_normal_alonge_wire_path.v(3,:),2);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% hold off;

end


end
