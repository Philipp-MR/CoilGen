function coil_parts= generate_cylindrical_pcb_print_(coil_parts,input)
%Generate a 2D pattern that can be rolled around a cylinder, @Philipp Amrein, Uniklinik
%Freiburg 2022
if input.surface_is_cylinder_flag & input.make_cylndrical_pcb

for part_ind=1:numel(coil_parts)

%Rotate the wire on the cylinder 
rot_mat=calc_3d_rotation_matrix_by_vector([input.cylinder_mesh_parameter_list(5) input.cylinder_mesh_parameter_list(6) input.cylinder_mesh_parameter_list(7)]',input.cylinder_mesh_parameter_list(8));
aligned_wire_path=(rot_mat*coil_parts.wire_path.v);

%for each point calculate the phi angle
phi_coord=atan2(aligned_wire_path(2,:),aligned_wire_path(1,:));

%open the wire on the end of angular wrap-arround; include cut-point for
%clean opening



end    

end

end
