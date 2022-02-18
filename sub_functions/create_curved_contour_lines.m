function contour_lines_curved=create_curved_contour_lines(contour_lines,planary_mesh,curved_mesh)
%Create the contour lines on the curved surface
contour_lines_curved=contour_lines;
%contour_lines_curved = rmfield(contour_lines_curved,'uvcoords');
% for loop_num=1:numel(contour_lines)
% contour_lines_curved(loop_num).point_coordinates=zeros(3,size(contour_lines(loop_num).uvcoords,2));
% end
for loop_num=1:numel(contour_lines)
contour_lines_curved(loop_num).point_coordinates=uv_to_xyz_loc(contour_lines(loop_num).uvcoords,planary_mesh,curved_mesh);
%disp(strcat('Loop Number',' ',num2str(loop_num),' Point Number ',num2str(point_num)));      
end

function points_out_3d=uv_to_xyz_loc(points_in_2d,planary_mesh,curved_mesh)
% from the 2D surface coordinates to the xyz coordinates of the 3D coil
% surface
points_out_3d=zeros(3,size(points_in_2d,2));
for point_ind=1:size(points_in_2d,2)
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,point_ind),points_in_2d(2,point_ind));
while isnan(target_triangle) % in case the point is directly between two triangles and no directly triangle indice can be found
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,point_ind)+ eps('double')*10*(1-2*rand(1)),points_in_2d(2,point_ind)+ eps('double')*10*(1-2*rand(1)));
end
points_out_3d(:,point_ind) = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
end
end

end
