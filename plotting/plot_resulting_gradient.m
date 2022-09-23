function  plot_resulting_gradient(coil_layouts,single_ind_to_plot,plot_title)

pos_data=coil_layouts(single_ind_to_plot).out.target_field.coords;

dot_size=100;

%Plot Histograms of the gradient
figure('name',plot_title);
tiledlayout('flow');
nexttile;
hold on
axis equal tight;
title("G"+" "+coil_layouts(single_ind_to_plot).out.input_data.field_shape_function+"[mT/m/A]",'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.gradient_field.*1000; % in mT\m\A
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
for part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts)
if isfield(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind),'wire_path')
plot3(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(2,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(3,:),'linewidth',2,'color',[0 0.4470 0.7410]);
else
for loop_ind=1:numel(coil_layouts.out.coil_parts(part_ind).contour_lines) 
plot3(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).contour_lines(loop_ind).v(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).contour_lines(loop_ind).v(2,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).contour_lines(loop_ind).v(3,:),'linewidth',2,'color',[0 0.4470 0.7410]);
end
end
end
caxis([coil_layouts(single_ind_to_plot).out.layout_gradient.mean_gradient-coil_layouts(single_ind_to_plot).out.layout_gradient.std_gradient coil_layouts(single_ind_to_plot).out.layout_gradient.mean_gradient+coil_layouts(single_ind_to_plot).out.layout_gradient.std_gradient]*1000);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
%Plot also the target gradient defined by the field shape
nexttile;
hold on;
title("Target Gradient", 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.gradient_field.*1000,'BinLimits',[0 10]);
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m/A]');
hold off;
set(gcf,'color','w');

end

