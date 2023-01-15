function planar_mesh=build_planar_mesh(planar_height,planar_width,num_lateral_divisions,num_longitudinal_divisions,rotation_vector_x,rotation_vector_y,rotation_vector_z,rotation_angle,center_position_x,center_position_y,center_position_z)
%create a mono-planar regular mesh in any orientation
% @Philipp Amrein, 2022, Uniklinik Freiburg

x_positions=-planar_width/2:(planar_width/(num_lateral_divisions)):planar_width/2;
y_positions=-planar_height/2:(planar_height/(num_longitudinal_divisions)):planar_height/2;

x=repmat(x_positions,[numel(y_positions) 1])';
y=repmat(y_positions',[1 numel(x_positions)])';

%define the vertic positions
vertices=calc_3d_rotation_matrix_by_vector([rotation_vector_x rotation_vector_y rotation_vector_z]',rotation_angle)*[y(:)'; x(:)'; zeros(size(y(:)'))]+[center_position_x center_position_y center_position_z]';

%define the mesh triangles
tri_1_vert_inds_1=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_1_vert_inds_2=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_1_vert_inds_3=[repmat(2:(num_lateral_divisions+1),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)-1]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_2_vert_inds_1=[tri_1_vert_inds_1];
tri_2_vert_inds_2=[repmat(1:(num_lateral_divisions),[num_longitudinal_divisions 1])+repmat([(1:num_longitudinal_divisions)]',[1 num_lateral_divisions])*(num_lateral_divisions+1)]';
tri_2_vert_inds_3=[tri_1_vert_inds_2];
faces_1=[tri_1_vert_inds_1(:) tri_1_vert_inds_2(:) tri_1_vert_inds_3(:)];
faces_2=[tri_2_vert_inds_1(:) tri_2_vert_inds_2(:) tri_2_vert_inds_3(:)];

planar_mesh.faces=[faces_1;faces_2];
planar_mesh.vertices=vertices';

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