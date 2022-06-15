function coil_parts= generate_cylindrical_pcb_print_(coil_parts,input)
%Generate a 2D pattern that can be rolled around a cylinder, @Philipp Amrein, Uniklinik
%Freiburg 2022

angular_shrink_factor=0.95;

if input.surface_is_cylinder_flag & input.make_cylndrical_pcb

if strcmp(input.pcb_interconnection_method,'spiral_in_out') % first layer spiral in; second layer spiral out


for part_ind=1:numel(coil_parts)

%Rotate the wire on the cylinder 
rot_mat=calc_3d_rotation_matrix_by_vector([input.cylinder_mesh_parameter_list(5) input.cylinder_mesh_parameter_list(6) input.cylinder_mesh_parameter_list(7)]',input.cylinder_mesh_parameter_list(8));
aligned_wire_path=(rot_mat*coil_parts.wire_path.v);

%for each point calculate the phi angle
phi_coord=atan2(aligned_wire_path(2,:),aligned_wire_path(1,:));

layout_2d=[phi_coord; aligned_wire_path(3,:)];

%open the wire on the end of angular wrap-arround;
cut_rectangle=[-pi*angular_shrink_factor pi*angular_shrink_factor pi*angular_shrink_factor -pi*angular_shrink_factor; ...
                                        max(aligned_wire_path(3,:))+abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        max(aligned_wire_path(3,:))+abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        min(aligned_wire_path(3,:))-abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        min(aligned_wire_path(3,:))-abs(max(aligned_wire_path(3,:))).*0.1];

cut_rectangle=[cut_rectangle cut_rectangle(:,1)];

%mark the segments of the wire which have angular wraparounds
wrap_edge=find(abs(diff(layout_2d(1,:)))>1.75*pi);

 %include cut-point for clean opening
for wrap_ind=1:numel(wrap_edge)

test_segment=[layout_2d(1,wrap_edge(wrap_ind)) layout_2d(1,wrap_edge(wrap_ind)+1); layout_2d(2,wrap_edge(wrap_ind)) layout_2d(2,wrap_edge(wrap_ind)+1)];
intersection_points=find_segment_intersections(test_segment,cut_rectangle);

%add the cut_points to the wire path

end

%rotate the alinged wire path around the z-axis in order to minimize the
%number of cuts which are neccessary around the angular wrap-arround


end    

else % first windings; second layer return paths


end

end

end
