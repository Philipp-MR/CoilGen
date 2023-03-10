function biplanar_mesh=build_biplanar_mesh(planar_height,planar_width,num_lateral_divisions,num_longitudinal_divisions,target_normal_x,target_normal_y,target_normal_z,center_position_x,center_position_y,center_position_z,plane_distance)
%create a mono-planar regular mesh in any orientation
% @Philipp Amrein, 2022, Uniklinik Freiburg


%define the mesh triangles
tri_1_vert_inds_1=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_1_vert_inds_2=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_1_vert_inds_3=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_2_vert_inds_1=[tri_1_vert_inds_1];
tri_2_vert_inds_2=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_2_vert_inds_3=[tri_1_vert_inds_2];
faces_1=[tri_1_vert_inds_1(:) tri_1_vert_inds_2(:) tri_1_vert_inds_3(:)];
faces_2=[tri_2_vert_inds_1(:) tri_2_vert_inds_2(:) tri_2_vert_inds_3(:)];

%define the vertex positions
x_positions=-planar_width/2:(planar_width/(num_lateral_divisions)):planar_width/2;
y_positions=-planar_height/2:(planar_height/(num_longitudinal_divisions)):planar_height/2;

x=repmat(x_positions,[numel(y_positions) 1])';
y=repmat(y_positions',[1 numel(x_positions)])';

old_normal=[0 0 1]';
target_normal=[target_normal_x,target_normal_y,target_normal_z]';
if norm(cross(old_normal,target_normal))~=0
rot_vec=cross(old_normal,target_normal)./norm(cross(old_normal,target_normal));
rot_angle = asin(norm(cross(old_normal,target_normal))/(norm(old_normal)*norm(target_normal)));
else
rot_vec=[1 0 0]';
rot_angle = 0;
end
 

z1=zeros(size(y))+plane_distance/2;
z2=zeros(size(y))-plane_distance/2;

rot_mat=calc_3d_rotation_matrix_by_vector([rot_vec(1) rot_vec(2) rot_vec(3)]',rot_angle);



vertices1=rot_mat*[y(:)'; x(:)'; z1(:)']+[center_position_x center_position_y center_position_z]';
vertices2=rot_mat*[y(:)'; x(:)'; z2(:)']+[center_position_x center_position_y center_position_z]';

faces_first_plane=[faces_1;faces_2];
faces_second_plane=[faces_1;faces_2]+size(vertices1,2);
% adjust the orientation of the second plane; to have symmetric normals
%faces_second_plane=[faces_second_plane(:,2) faces_second_plane(:,1) faces_second_plane(:,3)];

biplanar_mesh.faces=[faces_first_plane;faces_second_plane];
biplanar_mesh.vertices=[vertices1'; vertices2'];

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





%OLD version

% function biplanar_mesh=build_biplanar_mesh(planar_height,planar_width,num_lateral_divisions,num_longitudinal_divisions,rotation_vector_x,rotation_vector_y,rotation_vector_z,rotation_angle,center_position_x,center_position_y,center_position_z,plane_distance)
% %create a mono-planar regular mesh in any orientation
% % @Philipp Amrein, 2022, Uniklinik Freiburg
% 
% 
% %define the mesh triangles
% tri_1_vert_inds_1=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
% tri_1_vert_inds_2=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
% tri_1_vert_inds_3=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
% tri_2_vert_inds_1=[tri_1_vert_inds_1];
% tri_2_vert_inds_2=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
% tri_2_vert_inds_3=[tri_1_vert_inds_2];
% faces_1=[tri_1_vert_inds_1(:) tri_1_vert_inds_2(:) tri_1_vert_inds_3(:)];
% faces_2=[tri_2_vert_inds_1(:) tri_2_vert_inds_2(:) tri_2_vert_inds_3(:)];
% 
% %define the vertex positions
% x_positions=-planar_width/2:(planar_width/(num_lateral_divisions)):planar_width/2;
% y_positions=-planar_height/2:(planar_height/(num_longitudinal_divisions)):planar_height/2;
% 
% x=repmat(x_positions,[numel(y_positions) 1])';
% y=repmat(y_positions',[1 numel(x_positions)])';
% z1=zeros(size(y))+plane_distance/2;
% z2=zeros(size(y))-plane_distance/2;
% 
% rot_mat=calc_3d_rotation_matrix_by_vector([rotation_vector_x rotation_vector_y rotation_vector_z]',rotation_angle);
% vertices1=rot_mat*[y(:)'; x(:)'; z1(:)']+[center_position_x center_position_y center_position_z]';
% vertices2=rot_mat*[y(:)'; x(:)'; z2(:)']+[center_position_x center_position_y center_position_z]';
% 
% faces_first_plane=[faces_1;faces_2];
% faces_second_plane=[faces_1;faces_2]+size(vertices1,2);
% % adjust the orientation of the second plane; to have symmetric normals
% faces_second_plane=[faces_second_plane(:,2) faces_second_plane(:,1) faces_second_plane(:,3)];
% 
% biplanar_mesh.faces=[faces_first_plane;faces_second_plane];
% biplanar_mesh.vertices=[vertices1'; vertices2'];
% 
% function rot_mat_out= calc_3d_rotation_matrix_by_vector(rot_vec,rot_angle)
% %calculate the 3d rotation matrix around a rotation axis given by a vecor
% %and an angle
% %@Philipp Amrein 2022
% rot_vec=rot_vec./repmat(vecnorm(rot_vec),[3 1]); %normalize rot vector
% u_x=rot_vec(1);
% u_y=rot_vec(2);
% u_z=rot_vec(3);
% tmp1=sin(rot_angle);
% tmp2=cos(rot_angle);
% tmp3=(1-cos(rot_angle));
% rot_mat_out=zeros(3,3);
% rot_mat_out(1,1)=tmp2+u_x*u_x*tmp3;
% rot_mat_out(1,2)=u_x*u_y*tmp3-u_z*tmp1;
% rot_mat_out(1,3)=u_x*u_z*tmp3+u_y*tmp1;
% rot_mat_out(2,1)=u_y*u_x*tmp3+u_z*tmp1;
% rot_mat_out(2,2)=tmp2+u_y*u_y*tmp3;
% rot_mat_out(2,3)=u_y*u_z*tmp3-u_x*tmp1;
% rot_mat_out(3,1)=u_z*u_x*tmp3-u_y*tmp1;
% rot_mat_out(3,2)=u_z*u_y*tmp3+u_x*tmp1;
% rot_mat_out(3,3)=tmp2+u_z*u_z*tmp3;
% end
% 
% 
% end