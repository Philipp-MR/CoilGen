function  plot_various_error_metrics(coil_layouts,single_ind_to_plot,plot_title)



dot_size=200;
layout_c=coil_layouts(single_ind_to_plot).out.field_by_layout(3,:);
sf_c=coil_layouts(single_ind_to_plot).out.b_field_opt_sf(3,:);
loops_c=coil_layouts(single_ind_to_plot).out.field_by_unconnected_loops(3,:);
target_c=coil_layouts(single_ind_to_plot).out.target_field.b(3,:);
loops_c_1A=coil_layouts(single_ind_to_plot).out.field_loops_per1Amp(3,:);
layout_c_1A=coil_layouts(single_ind_to_plot).out.field_layout_per1Amp(3,:);
pos_data=coil_layouts(single_ind_to_plot).out.target_field.coords;

%Plot the target field in detail together with various error metrics

figure('name',plot_title);
tiledlayout('flow');
nexttile;
%Plot target field
hold on;
axis equal tight;
title('Target Bz, [mT/A]', 'interpreter', 'none');
plot_colors=target_c;
view(45,45);
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
h = colorbar;
ylabel(h, '[mT/A]');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;

%Plot field by stream function
hold on;
axis equal tight;
title('Bz by stream function, [mT/A]', 'interpreter', 'none');
plot_colors=sf_c;
view(45,45);
colormap(parula);
h = colorbar;
ylabel(h, '[mT/A]');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;
%Plot field by layout
hold on;
axis equal tight;
title('Layout Bz, [mT/A]', 'interpreter', 'none');
plot_colors=layout_c;
view(45,45);
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
h = colorbar;
ylabel(h, '[mT/A]');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;
%Plot field by unconnected contours
hold on;
axis equal tight;
title('Unconnected Contour Bz, [mT/A]', 'interpreter', 'none');
plot_colors=loops_c;
view(45,45);
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
h = colorbar;
ylabel(h, '[mT/A]');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off


nexttile;
%Plot relativ error between target and stream function field
hold on
axis equal tight;
title('Relative SF error, [%]', 'interpreter', 'none');
plot_colors=abs(sf_c-target_c)./max(abs(target_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;
%Plot relativ error between target and layout_c field
hold on
axis equal tight;
title('Relative error layout vs. target, [%]', 'interpreter', 'none');
plot_colors=abs(layout_c-target_c)./max(abs(target_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;
%Plot relativ error between layout and sf field
hold on
axis equal tight;
title('Relative error layout vs. sf field, [%]', 'interpreter', 'none');
plot_colors=abs(layout_c-sf_c)./max(abs(sf_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off


nexttile;
%Plot relativ error between target and unconnected contours
hold on
axis equal tight;
title('Relative error unconnected contours vs. target, [%]', 'interpreter', 'none');
plot_colors=abs(loops_c-target_c)./max(abs(target_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off

nexttile;
%Field differenz unconnected contours and final layout
hold on
axis equal tight;
title('Field difference between unconnected contours and final layout, [%]', 'interpreter', 'none');
plot_colors=abs(loops_c-layout_c)./max(abs(target_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off


set(gcf,'color','w');


end

