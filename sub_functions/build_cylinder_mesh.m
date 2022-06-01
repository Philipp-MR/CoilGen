function cylinder_mesh=build_cylinder_mesh(cylinder_height,cylinder_radius,num_circular_divisions,num_longitudinal_divisions,rotation_vector_x,rotation_vector_y,rotation_vector_z,rotation_angle)
%create a cylindrical regular mesh in any orientation
% @Philipp Amrein, 2022, Uniklinik Freiburg
x_positions=[sin(0:(2*pi)/(num_circular_divisions):2*pi)].*cylinder_radius;
y_positions=[cos(0:(2*pi)/(num_circular_divisions):2*pi)].*cylinder_radius;
x_positions(end)=[]; y_positions(end)=[]; %remove repetions at the ends
z_positions=((-1)*cylinder_height/2):(cylinder_height/(num_longitudinal_divisions+1)):(cylinder_height/2);

vertices_x=repmat(x_positions,[1 size(z_positions,2)-1]);
vertices_y=repmat(y_positions,[1 size(z_positions,2)-1]);
vertices_z=repelem(z_positions(1:end-1),size(x_positions,2));

vertices=[vertices_x; vertices_y; vertices_z];

tri_1_vert_inds_1=1:((num_circular_divisions)*(num_longitudinal_divisions));
tri_1_vert_inds_2=tri_1_vert_inds_1+1;
%take care of indice overflow at the end of the rings
tri_1_vert_inds_2(mod(tri_1_vert_inds_2-1,num_circular_divisions)==0)=tri_1_vert_inds_2(mod(tri_1_vert_inds_2-1,num_circular_divisions)==0)-num_circular_divisions;
tri_1_vert_inds_3=tri_1_vert_inds_2+num_circular_divisions;

tri_2_vert_inds_1=tri_1_vert_inds_1;
tri_2_vert_inds_2=tri_1_vert_inds_3;
tri_2_vert_inds_3=[1:((num_circular_divisions)*(num_longitudinal_divisions))]+num_circular_divisions;

faces_1=[tri_1_vert_inds_2; tri_1_vert_inds_1; tri_1_vert_inds_3];
faces_2=[tri_2_vert_inds_2; tri_2_vert_inds_1; tri_2_vert_inds_3];

cylinder_mesh.faces=[faces_1';faces_2']';
cylinder_mesh.vertices=vertices';


%rotate the cylinder in the desired orientation
rot_mat=calc_3d_rotation_matrix_by_vector([rotation_vector_x rotation_vector_y rotation_vector_z]',rotation_angle);
cylinder_mesh.vertices=(rot_mat*cylinder_mesh.vertices')';



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