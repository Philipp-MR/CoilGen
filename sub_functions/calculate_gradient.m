function layout_gradient= calculate_gradient(coil_parts,target_field,input)
%CALCULATE_LOCAL_GRADIENT
%calculate for each point the average of the local gradients to its neighbours


if strcmp(input.field_shape_function,'none')
field_function='z';
else
field_function=input.field_shape_function;
end

layout_gradient.dBxdxyz=zeros(size(target_field.coords));
layout_gradient.dBydxyz=zeros(size(target_field.coords));
layout_gradient.dBzdxyz=zeros(size(target_field.coords));
DBxdxyz=zeros(size(target_field.coords));
DBydxyz=zeros(size(target_field.coords));
DBzdxyz=zeros(size(target_field.coords));

if isfield(coil_parts(1),'wire_path')
for part_ind=1:numel(coil_parts)
[DBxdxyz,DBydxyz,DBzdxyz]=direct_biot_savart_gradient_calc(coil_parts(part_ind).wire_path.v,target_field.coords);
layout_gradient.dBxdxyz=layout_gradient.dBxdxyz+DBxdxyz;
layout_gradient.dBydxyz=layout_gradient.dBydxyz+DBydxyz;
layout_gradient.dBzdxyz=layout_gradient.dBzdxyz+DBzdxyz;
end
else %use the unconnected loops
for part_ind=1:numel(coil_parts)
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
[DBxdxyz,DBydxyz,DBzdxyz]=direct_biot_savart_gradient_calc(coil_parts(part_ind).contour_lines(loop_ind).v,target_field.coords);
layout_gradient.dBxdxyz=layout_gradient.dBxdxyz+DBxdxyz;
layout_gradient.dBydxyz=layout_gradient.dBydxyz+DBydxyz;
layout_gradient.dBzdxyz=layout_gradient.dBzdxyz+DBzdxyz;
end
end
end


my_fun=str2func("@(x,y,z)"+field_function);
norm_dir_x=my_fun(1,0,0);
norm_dir_y=my_fun(0,1,0);
norm_dir_z=my_fun(0,0,1);
gradient_direction=[norm_dir_x norm_dir_y norm_dir_z];
gradient_direction=gradient_direction./norm(gradient_direction);

%Project the gradient direction to the full set of cartesian gradients
layout_gradient.gradient_in_target_direction=((gradient_direction(1).*layout_gradient.dBzdxyz(1,:)).^2+(gradient_direction(2).*layout_gradient.dBzdxyz(2,:)).^2+(gradient_direction(3).*layout_gradient.dBzdxyz(3,:)).^2).^(1/2);
layout_gradient.mean_gradient_in_target_direction=mean(layout_gradient.gradient_in_target_direction,'omitnan');
layout_gradient.std_gradient_in_target_direction=std(layout_gradient.gradient_in_target_direction,'omitnan');

%Adjust to the unit [mT\m\A]
layout_gradient.dBxdxyz=layout_gradient.dBxdxyz.*1000;
layout_gradient.dBydxyz=layout_gradient.dBydxyz.*1000;
layout_gradient.dBzdxyz=layout_gradient.dBzdxyz.*1000;
layout_gradient.gradient_in_target_direction=layout_gradient.gradient_in_target_direction.*1000;
layout_gradient.mean_gradient_in_target_direction=layout_gradient.mean_gradient_in_target_direction.*1000;
layout_gradient.std_gradient_in_target_direction=layout_gradient.std_gradient_in_target_direction.*1000;


function [DBxdxyz,DBydxyz,DBzdxyz]=direct_biot_savart_gradient_calc(wire_elements,target_coords)
%Calculate b field with biot savarts law  by wire elements given as a sequence of coordinate points
num_tp=size(target_coords,2);
%to avoid memory problems splitt the curve into several parts
track_part_length=1000;
if size(wire_elements,2)>track_part_length
track_part_inds=[1:track_part_length:size(wire_elements,2) size(wire_elements,2)];
if track_part_inds(end-1)==track_part_inds(end)
track_part_inds(end)=[];
end
for parts_ind=1:numel(track_part_inds)-1
wire_part(parts_ind).coord=wire_elements(:,track_part_inds(parts_ind):track_part_inds(parts_ind+1));
wire_part(parts_ind).seg_coords=(wire_part(parts_ind).coord(:,1:end-1)+wire_part(parts_ind).coord(:,2:end))./2;
wire_part(parts_ind).currents=wire_part(parts_ind).coord(:,2:end)-wire_part(parts_ind).coord(:,1:end-1);
%wire_part(part_ind).seg_length=wire_elements(:,track_part_inds(part_ind):track_part_inds(part_ind+1));
end
else
wire_part.coord=wire_elements;
wire_part.seg_coords=(wire_part.coord(:,1:end-1)+wire_part.coord(:,2:end))./2;
wire_part.currents=wire_part.coord(:,2:end)-wire_part.coord(:,1:end-1);
end
%Calculate the magnetic field with Biot-Savarts law
DBxdxyz=zeros(3,num_tp);
DBydxyz=zeros(3,num_tp);
DBzdxyz=zeros(3,num_tp);
for parts_ind=1:numel(wire_part)
target_p=repmat(target_coords,[1 1 size(wire_part(parts_ind).seg_coords,2)]); % copies target point coordinates
target_p=permute(target_p,[1 3 2]);
cur_pos= repmat(wire_part(parts_ind).seg_coords(:,:),[1 1 num_tp]); % copies vector of current position
cur_dir= repmat(wire_part(parts_ind).currents(:,:),[1 1 num_tp]); % copies vector of current direction

%Calculate the gradient:
x=squeeze(target_p(1,:,:));
y=squeeze(target_p(2,:,:));
z=squeeze(target_p(3,:,:));
l_x=squeeze(cur_pos(1,:,:));
l_y=squeeze(cur_pos(2,:,:));
l_z=squeeze(cur_pos(3,:,:));
d_l_x=squeeze(cur_dir(1,:,:));
d_l_y=squeeze(cur_dir(2,:,:));
d_l_z=squeeze(cur_dir(3,:,:));

E=((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(3/2);
dEdx=3.*(x-l_x).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2);
dEdy=3.*(y-l_y).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2);
dEdz=3.*(z-l_z).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2); 
C=10^(-7);

phi_x=(d_l_y.*(z-l_z)-d_l_z.*(y-l_y));
phi_y=(d_l_z.*(x-l_x)-d_l_x.*(z-l_z));
phi_z=(d_l_x.*(y-l_y)-d_l_y.*(x-l_x));

theta_factor=((E.*E).^(-1)).*(-1).*C;

%gradient of bx compnent
dBx_dx=theta_factor.*dEdx.*phi_x;
%gradient in y of bz compnent
dBx_dy=theta_factor.*dEdy.*phi_x+(E.^(-1)).*C.*(-1).*d_l_z;
%gradient in z of bz compnent
dBx_dz=theta_factor.*dEdz.*phi_x+(E.^(-1)).*C.*d_l_y;

%gradient of by compnent
dBy_dx=theta_factor.*dEdx.*phi_y+(E.^(-1)).*C.*d_l_z;
%gradient in y of bz compnent
dBy_dy=theta_factor.*dEdy.*phi_y;
%gradient in z of bz compnent
dBy_dz=theta_factor.*dEdz.*phi_y+(E.^(-1)).*C.*(-1).*d_l_x;

%gradient of bz compnent
dBz_dx=theta_factor.*dEdx.*phi_z+(E.^(-1)).*C.*(-1).*d_l_y;
%gradient in y of bz compnent
dBz_dy=theta_factor.*dEdy.*phi_z+(E.^(-1)).*C.*d_l_x;
%gradient in z of bz compnent
dBz_dz=theta_factor.*dEdz.*phi_z;

DBxdxyz=DBxdxyz+[sum(dBx_dx,1); sum(dBx_dy,1); sum(dBx_dz,1)];
DBydxyz=DBydxyz+[sum(dBy_dx,1); sum(dBy_dy,1); sum(dBy_dz,1)];
DBzdxyz=DBzdxyz+[sum(dBz_dx,1); sum(dBz_dy,1); sum(dBz_dz,1)];
end
end

%OLD VERSION

% mesh_diag_length=norm([max(target_field.coords(1,:))-min(target_field.coords(1,:));max(target_field.coords(1,:))-min(target_field.coords(1,:));max(target_field.coords(1,:))-min(target_field.coords(1,:))]);
% gradient_out=calc_gradient_along_vector(coil_parts,target_field.coords,field_function,mesh_diag_length/10000);

% % function gradient_out=calc_gradient_along_vector(coil_parts,field_coords,target_endcoding_function,delta_shift_length)
% % %calculate the mean gradient in a given direction
% % %and an angle
% % %two versions of the target fields which are slightly shifted in their
% % %coordiantes in the direction of the aimed gradient
% % %the diffenerce of those fields allows than to obtain the local gradient
% % %for each point
% % %@Philipp Amrein 2022
% % 
% % my_fun=str2func("@(x,y,z)"+target_endcoding_function);
% % norm_dir_x=my_fun(1,0,0);
% % norm_dir_y=my_fun(0,1,0);
% % norm_dir_z=my_fun(0,0,1);
% % target_direction=[0 0 1];
% % gradient_direction=[norm_dir_x norm_dir_y norm_dir_z];
% % gradient_direction=gradient_direction./norm(gradient_direction);
% % if norm(cross(gradient_direction,target_direction))~=0
% % rot_vector=cross(gradient_direction,target_direction)./norm(cross(gradient_direction,target_direction));
% % rot_angle=asin(norm(cross(gradient_direction,target_direction))/(norm(target_direction)*norm(gradient_direction)));
% % else
% % rot_vector=[1 0 0];
% % rot_angle = 0;
% % end
% % rot_mat_out= local_calc_3d_rotation_matrix(rot_vector',rot_angle);
% % rotated_field_coords=rot_mat_out*(field_coords-mean(field_coords,2))+mean(field_coords,2);
% % %rotated_field_coords=rot_mat_out*field_coords;
% % delta_shift=[zeros(1,size(field_coords,2)); zeros(1,size(field_coords,2)); ones(1,size(field_coords,2))].*delta_shift_length;
% % %create two target fields slightly shifted in gradient direction 
% % target_fied_1=inv(rot_mat_out)*(rotated_field_coords-delta_shift);
% % target_fied_2=inv(rot_mat_out)*(rotated_field_coords+delta_shift);
% % for part_ind=1:numel(coil_parts)
% % if isfield(coil_parts(part_ind),'wire_path')
% % b1=local_biot_savart_calc_b(coil_parts(part_ind).wire_path.v,target_fied_1);
% % b2=local_biot_savart_calc_b(coil_parts(part_ind).wire_path.v,target_fied_2);
% % else %use the loops
% % b1=zeros(size(target_fied_1));
% % b2=zeros(size(target_fied_2));
% % for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
% % b1=b1+local_biot_savart_calc_b(coil_parts(part_ind).contour_lines(loop_ind).v,target_fied_1);
% % b2=b2+local_biot_savart_calc_b(coil_parts(part_ind).contour_lines(loop_ind).v,target_fied_2);
% % end
% % end
% % gradient_out=(b2(3,:)-b1(3,:))./delta_shift_length;
% % end
% % 
% % end
% % 
% % 
% % function rot_mat_out= local_calc_3d_rotation_matrix(rot_vec,rot_angle)
% % %calculate the 3d rotation matrix around a rotation axis given by a vecor
% % %and an angle
% % %@Philipp Amrein 2022
% % rot_vec=rot_vec./repmat(vecnorm(rot_vec),[3 1]); %normalize rot vector
% % u_x=rot_vec(1);
% % u_y=rot_vec(2);
% % u_z=rot_vec(3);
% % tmp1=sin(rot_angle);
% % tmp2=cos(rot_angle);
% % tmp3=(1-cos(rot_angle));
% % rot_mat_out=zeros(3,3);
% % rot_mat_out(1,1)=tmp2+u_x*u_x*tmp3;
% % rot_mat_out(1,2)=u_x*u_y*tmp3-u_z*tmp1;
% % rot_mat_out(1,3)=u_x*u_z*tmp3+u_y*tmp1;
% % rot_mat_out(2,1)=u_y*u_x*tmp3+u_z*tmp1;
% % rot_mat_out(2,2)=tmp2+u_y*u_y*tmp3;
% % rot_mat_out(2,3)=u_y*u_z*tmp3-u_x*tmp1;
% % rot_mat_out(3,1)=u_z*u_x*tmp3-u_y*tmp1;
% % rot_mat_out(3,2)=u_z*u_y*tmp3+u_x*tmp1;
% % rot_mat_out(3,3)=tmp2+u_z*u_z*tmp3;
% % end
% % 
% % function b_field=local_biot_savart_calc_b(wire_elements,target_coords)
% % %Calculate b field with biot savarts law  by wire elements given as a sequence of coordinate points
% % num_tp=size(target_coords,2);
% % %to avoid memory problems splitt the curve into several parts
% % track_part_length=1000;
% % if size(wire_elements,2)>track_part_length
% % track_part_inds=[1:track_part_length:size(wire_elements,2) size(wire_elements,2)];
% % if track_part_inds(end-1)==track_part_inds(end)
% % track_part_inds(end)=[];
% % end
% % for parts_ind=1:numel(track_part_inds)-1
% % wire_part(parts_ind).coord=wire_elements(:,track_part_inds(parts_ind):track_part_inds(parts_ind+1));
% % wire_part(parts_ind).seg_coords=(wire_part(parts_ind).coord(:,1:end-1)+wire_part(parts_ind).coord(:,2:end))./2;
% % wire_part(parts_ind).currents=wire_part(parts_ind).coord(:,2:end)-wire_part(parts_ind).coord(:,1:end-1);
% % %wire_part(part_ind).seg_length=wire_elements(:,track_part_inds(part_ind):track_part_inds(part_ind+1));
% % end
% % else
% % wire_part.coord=wire_elements;
% % wire_part.seg_coords=(wire_part.coord(:,1:end-1)+wire_part.coord(:,2:end))./2;
% % wire_part.currents=wire_part.coord(:,2:end)-wire_part.coord(:,1:end-1);
% % end
% % %Calculate the magnetic field with Biot-Savarts law
% % b_field=zeros(3,num_tp);
% % for parts_ind=1:numel(wire_part)
% % target_p=repmat(target_coords,[1 1 size(wire_part(parts_ind).seg_coords,2)]); % copies target point coordinates
% % target_p=permute(target_p,[1 3 2]);
% % cur_pos= repmat(wire_part(parts_ind).seg_coords(:,:),[1 1 num_tp]); % copies vector of current position
% % cur_dir= repmat(wire_part(parts_ind).currents(:,:),[1 1 num_tp]); % copies vector of current direction
% % R = target_p-cur_pos;  % distance from current path to each point
% % %R(:,sqrt(dot(R,R,1))<min_distance_to_wire)=10^12; % if distance is to small the distance is blown up to ignore large values
% % len_R=(1./vecnorm(R,2,1)).^3; % distance factor of biot savart law
% % len_R=repmat(len_R,[3 1 1]);
% % dB=10^(-7)*cross(cur_dir,R,1).*len_R;
% % b_field(:,:)=b_field(:,:)+squeeze(sum(dB,2));
% % end
% % end
% % 
% % %Another idea would be to select the unique x,y,z values of the point cloud
% % %and build a meshgrid on it for which the gradient can then be regularly
% % %calculated ...

end

