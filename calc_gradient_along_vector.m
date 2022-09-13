function [mean_gradient_strength,gradient_out]=calc_gradient_along_vector(field,field_coords,target_endcoding_function)
%calculate the mean gradient in a given direction
%and an angle
%@Philipp Amrein 2022

my_fun=str2func("@(x,y,z)"+target_endcoding_function);

norm_dir_x=my_fun(1,0,0);
norm_dir_y=my_fun(0,1,0);
norm_dir_z=my_fun(0,0,1);

target_direction=[0 0 1];
gradient_direction=[norm_dir_x norm_dir_y norm_dir_z];
gradient_direction=gradient_direction./norm(gradient_direction);

if norm(cross(gradient_direction,target_direction))~=0
rot_vector=cross(gradient_direction,target_direction)./norm(cross(gradient_direction,target_direction));
rot_angle=asin(norm(cross(gradient_direction,target_direction))/(norm(target_direction)*norm(gradient_direction)));
else
rot_vector=[1 0 0];
rot_angle = 0;
end



rot_mat_out= calc_3d_rotation_matrix(rot_vector',rot_angle);
rotated_field_coords=rot_mat_out*(field_coords-mean(field_coords,2));
gradient_out=field(3,:)./rotated_field_coords(3,:);
gradient_out(abs(rotated_field_coords(3,:))<10^(-6))=nan;
mean_gradient_strength=mean(gradient_out,'omitnan');


function rot_mat_out= calc_3d_rotation_matrix(rot_vec,rot_angle)
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