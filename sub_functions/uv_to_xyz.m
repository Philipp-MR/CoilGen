function [points_out_3d,points_in_2d]=uv_to_xyz(points_in_2d,planary_mesh,curved_mesh)
% from the 2D surface coordinates to the xyz coordinates of the 3D coil
% surface
num_deleted_points=0;
points_out_3d=zeros(3,size(points_in_2d,2));
avg_mesh_diameter=norm([planary_mesh.Points-mean(planary_mesh.Points,1)]');


for rrrr=1:size(points_in_2d,2)
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,rrrr-num_deleted_points),points_in_2d(2,rrrr-num_deleted_points));
tri_inds=0;
while isnan(target_triangle) % in case the point is directly between two triangles and no directly triangle indice can be found
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,rrrr-num_deleted_points)+ avg_mesh_diameter*(0.5-rand(1))/1000,points_in_2d(2,rrrr-num_deleted_points)+avg_mesh_diameter*(0.5-rand(1))/1000);
tri_inds=tri_inds+1;
if tri_inds>1000
    fprintf(' Warning: Points cannot be assigned to a triangle \n');
    break;
end
end
if ~isnan(target_triangle)
points_out_3d(:,rrrr-num_deleted_points) = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
else %remove the point
points_in_2d(:,rrrr-num_deleted_points)=[];
points_out_3d(:,rrrr-num_deleted_points)=[];
num_deleted_points=num_deleted_points+1;
end
end
end