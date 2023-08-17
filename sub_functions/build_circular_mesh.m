function circular_mesh=build_circular_mesh(radius,num_radial_divisions,rotation_vector_x,rotation_vector_y,rotation_vector_z,rotation_angle,center_position_x,center_position_y,center_position_z)
%create a mono-planar circular mesh in any orientation
% @Philipp Amrein, 2022, Uniklinik Freiburg

x_positions=-radius:(radius*2/(num_radial_divisions)):radius;
y_positions=-radius:(radius*2/(num_radial_divisions)):radius;

x=repmat(x_positions,[numel(y_positions) 1])';
y=repmat(y_positions',[1 numel(x_positions)])';

%define the vertic positions
vertices=[y(:)'; x(:)'; zeros(size(y(:)'))];


%define the mesh triangles
tri_1_vert_inds_1=[repmat(1:(num_radial_divisions),[num_radial_divisions 1])+repmat([(1:num_radial_divisions)-1]',[1 num_radial_divisions])*(num_radial_divisions+1)]';
tri_1_vert_inds_2=[repmat(2:(num_radial_divisions+1),[num_radial_divisions 1])+repmat([(1:num_radial_divisions)]',[1 num_radial_divisions])*(num_radial_divisions+1)]';
tri_1_vert_inds_3=[repmat(2:(num_radial_divisions+1),[num_radial_divisions 1])+repmat([(1:num_radial_divisions)-1]',[1 num_radial_divisions])*(num_radial_divisions+1)]';
tri_2_vert_inds_1=[tri_1_vert_inds_1];
tri_2_vert_inds_2=[repmat(1:(num_radial_divisions),[num_radial_divisions 1])+repmat([(1:num_radial_divisions)]',[1 num_radial_divisions])*(num_radial_divisions+1)]';
tri_2_vert_inds_3=[tri_1_vert_inds_2];
faces_1=[tri_1_vert_inds_1(:) tri_1_vert_inds_2(:) tri_1_vert_inds_3(:)];
faces_2=[tri_2_vert_inds_1(:) tri_2_vert_inds_2(:) tri_2_vert_inds_3(:)];

circular_mesh.faces=[faces_1;faces_2];
circular_mesh.vertices=vertices';


%Delete ther vertices outside the ciruclar boundary
vertices_to_delete=sort(find((circular_mesh.vertices(:,1).^2+circular_mesh.vertices(:,2).^2+circular_mesh.vertices(:,3).^2).^(1/2)>radius*0.99),'descend');
for delete_ind=1:numel(vertices_to_delete)
circular_mesh.vertices(vertices_to_delete(delete_ind),:)=[];
%remove faces related to the deleted vertex and update the index list
faces_to_delete=find(circular_mesh.faces(:,1)==vertices_to_delete(delete_ind)|circular_mesh.faces(:,2)==vertices_to_delete(delete_ind)|circular_mesh.faces(:,3)==vertices_to_delete(delete_ind));
circular_mesh.faces(faces_to_delete,:)=[];
circular_mesh.faces(circular_mesh.faces>vertices_to_delete(delete_ind))=circular_mesh.faces(circular_mesh.faces>vertices_to_delete(delete_ind))-1;
end

% Morph the boundary verts for proper circular shape
boundary_verts=freeBoundary(triangulation(circular_mesh.faces,circular_mesh.vertices));
boundary_verts=unique(boundary_verts(:));
for vert_ind=1:numel(boundary_verts)
vert_radius=(circular_mesh.vertices(boundary_verts(vert_ind),1).^2+circular_mesh.vertices(boundary_verts(vert_ind),2).^2).^(1/2);
circular_mesh.vertices(boundary_verts(vert_ind),:)=circular_mesh.vertices(boundary_verts(vert_ind),:).*(radius/vert_radius);
end


% %Add addiational faces&vertices for proper circular boundary
% num_end_ring_vertices=ceil(num_radial_divisions*pi);
% end_ring_x=sin(0:(2*pi)/(num_end_ring_vertices-1):(2*pi)).*radius;
% end_ring_y=cos(0:(2*pi)/(num_end_ring_vertices-1):(2*pi)).*radius;
% end_ring_z=zeros(size(end_ring_x));
% end_ring_vertices=[end_ring_x(1:end-1); end_ring_y(1:end-1); end_ring_z(1:end-1)];
% first_verts=1:(numel(end_ring_vertices)-1);
% second_verts=2:numel(end_ring_vertices);
% used_verts=zeros(1,size(planar_mesh.vertices,1));
% dists=[];
% for vert_ind=1:numel(end_ring_vertices)-1
% dists=((planar_mesh.vertices(:,1)-end_ring_vertices(1,vert_ind)).^2+(planar_mesh.vertices(:,2)-end_ring_vertices(2,vert_ind)).^2+(planar_mesh.vertices(:,3)-end_ring_vertices(3,vert_ind)).^2).^(1/2);
% [~,sort_dist_inds] = sort(dists);
% face1_vert_1=first_verts(vert_ind);
% face1_vert_2=second_verts(vert_ind);
% face1_vert_3=sort_dist_inds(1);
% face_ind=vert_ind+1;
% end

%Adjust the mesh in the with desired translation and rotation
circular_mesh.vertices=[calc_3d_rotation_matrix_by_vector([rotation_vector_x rotation_vector_y rotation_vector_z]',rotation_angle)*circular_mesh.vertices'+[center_position_x center_position_y center_position_z]']';

% figure; 
% hold on;
% trisurf(triangulation(planar_mesh.faces,planar_mesh.vertices),'facecolor','cyan');
% scatter3(end_ring_vertices(1,:),end_ring_vertices(2,:),end_ring_vertices(3,:));
% plot3(end_ring_x,end_ring_y,end_ring_z,'r','linewidth',2);
% %axis equal; 
% view(45,45)
% hold off;




function rot_mat_out= calc_3d_rotation_matrix_by_vector(rot_vec,rot_angle)
%calculate the 3d rotation matrix around a rotation axis given by a vecor
%and an angle
%@Philipp Amrein 2022
rot_vec=rot_vec./repmat(vecnorm(rot_vec),[3 1]); %normalize rot vector
u_x=rot_vec(1);
u_y=rot_vec(2);
u_z=rot_vec(3);
tmp1=sin(rot_angle);
tmp2=cos(rot_angle);
tmp3=(1-cos(rot_angle));
rot_mat_out=zeros(3,3);
rot_mat_out(1,1)=tmp2+u_x*u_x*tmp3;
rot_mat_out(1,2)=u_x*u_y*tmp3-u_z*tmp1;
rot_mat_out(1,3)=u_x*u_z*tmp3+u_y*tmp1;
rot_mat_out(2,1)=u_y*u_x*tmp3+u_z*tmp1;
rot_mat_out(2,2)=tmp2+u_y*u_y*tmp3;
rot_mat_out(2,3)=u_y*u_z*tmp3-u_x*tmp1;
rot_mat_out(3,1)=u_z*u_x*tmp3-u_y*tmp1;
rot_mat_out(3,2)=u_z*u_y*tmp3+u_x*tmp1;
rot_mat_out(3,3)=tmp2+u_z*u_z*tmp3;
end



end