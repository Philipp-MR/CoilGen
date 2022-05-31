function local_opening_gab=calc_local_opening_gab(loop,point_1,point_2,opening_gab)

%local_opening_gab=calc_local_opening_gab2(coil_mesh,loop,cut_point,cut_direction,opening_gab)

if ~isempty(point_2) %two points are specified

uv_distance=vecnorm(loop.uv(:,point_1)-loop.uv(:,point_2));
v_distane=vecnorm(loop.v(:,point_1)-loop.v(:,point_2));

local_opening_gab=opening_gab*uv_distance/v_distane;

else %only one point is specified, find the other to build the direction

[~,min_ind_2]=min(vecnorm(loop.uv-point_1));

min_ind_1=min_ind_2+2;

if min_ind_1<0; min_ind_1=min_ind_1+size(loop.uv,2); end;
if min_ind_1>size(loop.uv,2); min_ind_1=min_ind_1-size(loop.uv,2); end;

uv_distance=vecnorm(loop.uv(:,min_ind_1)-loop.uv(:,min_ind_2));
v_distane=vecnorm(loop.v(:,min_ind_1)-loop.v(:,min_ind_2));

local_opening_gab=opening_gab*uv_distance/v_distane;


end


%Old version with edge projection
%Calulate the local opening width within the flat 2D domian
% target_triangle=pointLocation(triangulation(coil_mesh.f,coil_mesh.uv'),cut_point.uv(1),cut_point.uv(2));
% 
% vertex(1).v=coil_mesh.v(coil_mesh.f(target_triangle,1),:)';
% vertex(1).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,1));
% 
% vertex(2).v=coil_mesh.v(coil_mesh.f(target_triangle,2),:)';
% vertex(2).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,2));
% 
% vertex(3).v=coil_mesh.v(coil_mesh.f(target_triangle,3),:)';
% vertex(3).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,3));
% 
% edge(1).v=vertex(2).v-vertex(1).v;
% edge(2).v=vertex(3).v-vertex(2).v;
% edge(3).v=vertex(1).v-vertex(3).v;
% 
% edge(1).uv=vertex(2).uv-vertex(1).uv;
% edge(2).uv=vertex(3).uv-vertex(2).uv;
% edge(3).uv=vertex(1).uv-vertex(3).uv;
% 
% %calculate the length of shift in the wanted direction
% cut_direction=cut_direction./vecnorm(cut_direction);
% 
% cut_direction_xyz=dot(cut_direction,edge(1).uv)*edge(1).v+dot(cut_direction,edge(2).uv)*edge(2).v;
% curved_unit_length=vecnorm(cut_direction_xyz);
% 
% % calcalte the wanted length that corresponds to the targeted opening widht
% % but in the uv domain
% local_opening_gab=opening_gab/curved_unit_length;








% % %Old version with Transformtion Matrix
% % 
% % %Calulate the local opening width within the flat 2D domian
% % target_triangle=pointLocation(triangulation(coil_mesh.f,coil_mesh.uv'),cut_point.uv(1),cut_point.uv(2));
% % 
% % vertex(1).v=coil_mesh.v(coil_mesh.f(target_triangle,1),:)';
% % vertex(1).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,1));
% % 
% % vertex(2).v=coil_mesh.v(coil_mesh.f(target_triangle,2),:)';
% % vertex(2).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,2));
% % 
% % vertex(3).v=coil_mesh.v(coil_mesh.f(target_triangle,3),:)';
% % vertex(3).uv=coil_mesh.uv(:,coil_mesh.f(target_triangle,3));
% % 
% % 
% % % all_dx=[vertex(2).v(1)-vertex(1).v(1) vertex(3).v(1)-vertex(1).v(1) vertex(3).v(1)-vertex(2).v(1)];
% % % all_dy=[vertex(2).v(2)-vertex(1).v(2) vertex(3).v(2)-vertex(1).v(2) vertex(3).v(2)-vertex(2).v(2)];
% % % all_dz=[vertex(2).v(3)-vertex(1).v(3) vertex(3).v(3)-vertex(1).v(3) vertex(3).v(3)-vertex(2).v(3)];
% % 
% % all_du=[vertex(2).uv(1)-vertex(1).uv(1) vertex(3).uv(1)-vertex(1).uv(1) vertex(3).v(1)-vertex(2).uv(1)];
% % all_dv=[vertex(2).uv(2)-vertex(1).uv(2) vertex(3).uv(2)-vertex(1).uv(2) vertex(3).v(2)-vertex(2).uv(2)];
% % 
% % %find the vertices with largest coordinate differnence
% % 
% % % [~,dx_max_ind]=max((all_dx));
% % % [~,dy_max_ind]=max((all_dy));
% % % [~,dz_max_ind]=max((all_dz));
% % % 
% % % [~,dx_min_ind]=min((all_dx));
% % % [~,dy_min_ind]=min((all_dy));
% % % [~,dz_min_ind]=min((all_dz));
% % 
% % [~,du_max_ind]=max((all_du));
% % [~,dv_max_ind]=max((all_dv));
% % 
% % [~,du_min_ind]=min((all_du));
% % [~,dv_min_ind]=min((all_dv));
% % 
% % 
% % % % %Calculate the elements of the locally linear xyz->uv transformation matrix
% % % du_dx=(vertex(dx_max_ind).uv(1)-vertex(dx_min_ind).uv(1))/(vertex(dx_max_ind).v(1)-vertex(dx_min_ind).v(1));
% % % du_dy=(vertex(dy_max_ind).uv(1)-vertex(dy_min_ind).uv(1))/(vertex(dy_max_ind).v(2)-vertex(dy_min_ind).v(2));
% % % du_dz=(vertex(dz_max_ind).uv(1)-vertex(dz_min_ind).uv(1))/(vertex(dz_max_ind).v(3)-vertex(dz_min_ind).v(3));
% % % 
% % % dv_dx=(vertex(dx_max_ind).uv(2)-vertex(dx_min_ind).uv(2))/(vertex(dx_max_ind).v(1)-vertex(dx_min_ind).v(1));
% % % dv_dy=(vertex(dy_max_ind).uv(2)-vertex(dy_min_ind).uv(2))/(vertex(dy_max_ind).v(2)-vertex(dy_min_ind).v(2));
% % % dv_dz=(vertex(dz_max_ind).uv(2)-vertex(dz_min_ind).uv(2))/(vertex(dz_max_ind).v(3)-vertex(dz_min_ind).v(3));
% % 
% % %Calculate the elements of the locally linear uv->xyz transformation matrix
% % dx_du=(vertex(du_max_ind).v(1)-vertex(du_min_ind).v(1))/(vertex(du_max_ind).uv(1)-vertex(du_min_ind).uv(1));
% % dx_dv=(vertex(dv_max_ind).v(1)-vertex(dv_min_ind).v(1))/(vertex(dv_max_ind).uv(2)-vertex(dv_min_ind).uv(2));
% % 
% % dy_du=(vertex(du_max_ind).v(2)-vertex(du_min_ind).v(2))/(vertex(du_max_ind).uv(1)-vertex(du_min_ind).uv(1));
% % dy_dv=(vertex(dv_max_ind).v(2)-vertex(dv_min_ind).v(2))/(vertex(dv_max_ind).uv(2)-vertex(dv_min_ind).uv(2));
% % 
% % dz_du=(vertex(du_max_ind).v(3)-vertex(du_min_ind).v(3))/(vertex(du_max_ind).uv(1)-vertex(du_min_ind).uv(1));
% % dz_dv=(vertex(dv_max_ind).v(3)-vertex(dv_min_ind).v(3))/(vertex(dv_max_ind).uv(2)-vertex(dv_min_ind).uv(2));
% % 
% % %Build the transformation matrix
% % uv_to_xyz_mat=[dx_du dx_dv; dy_du dy_dv; dz_du dz_dv];
% % % xyz_to_uv_mat=[du_dx du_dy du_dz; dv_dx dv_dy dv_dz];
% % 
% % 
% % %calculate the length of shift in the wanted direction
% % cut_direction=cut_direction./vecnorm(cut_direction);
% % shift_vec=uv_to_xyz_mat*cut_direction;
% % curved_unit_length=vecnorm(shift_vec);
% % 
% % %calcalte the wanted length that corresponds to the targeted opening widht
% % %but in the uv domain
% % local_opening_gab=opening_gab/curved_unit_length;





end