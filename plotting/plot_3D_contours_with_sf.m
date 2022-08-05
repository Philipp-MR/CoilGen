function  plot_3D_contours_with_sf(coil_layout,single_ind_to_plot,plot_title)

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
    fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',[normed_pot_1 normed_pot_2 normed_pot_3]');
    if isfield(coil_layout(single_ind_to_plot).out.coil_parts(1),'groups')
    for group_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups)
    for contour_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops)
    plot3(  coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(contour_ind).v(1,:)...
                ,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(contour_ind).v(2,:)...
                ,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(contour_ind).v(3,:),'color','k','linewidth',2);
    end
    end
    end
end
% %Together with target field
% plot_colors=coil_layout(single_ind_to_plot).out.target_field.b(3,:)./max(abs(coil_layout(single_ind_to_plot).out.target_field.b(3,:)));
% scatter3(coil_layout(single_ind_to_plot).out.target_field.coords(1,:),...
%     coil_layout(single_ind_to_plot).out.target_field.coords(2,:),...
%     coil_layout(single_ind_to_plot).out.target_field.coords(3,:),...
%     100*ones(1,numel(coil_layout(single_ind_to_plot).out.target_field.coords,2)),...
%     plot_colors,'filled');
% axis equal tight;
% caxis([min(coil_layout.out.combined_mesh.stream_function) max(coil_layout.out.combined_mesh.stream_function)]);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
set(gcf,'color','w');
axis equal;
hold off;


end

