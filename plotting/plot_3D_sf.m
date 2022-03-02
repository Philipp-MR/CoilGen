function  plot_3D_sf(coil_layout,single_ind_to_plot,plot_title)

% Plot the stream function interpolatet on triangular mesh

figure;
hold on;
title(plot_title+": "+'Stream function by optimization and target Bz', 'interpreter', 'none');
view(3);
for part_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts)
    normed_pot=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).stream_function;
    normed_pot=normed_pot./max(abs(normed_pot));
	x_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    x_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    x_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    y_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    y_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    y_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    z_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    z_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    z_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    normed_pot_1=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    normed_pot_2=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    normed_pot_3=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',[normed_pot_1 normed_pot_2 normed_pot_3]','EdgeAlpha',0);
end

axis equal;
hold off;

end

