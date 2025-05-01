function  plot_various_error_metrics(coil_layouts,single_ind_to_plot,plot_title)



dot_size=200;
layout_c=coil_layouts(single_ind_to_plot).out.field_layout_per1Amp(3,:).*1000;
sf_c=coil_layouts(single_ind_to_plot).out.b_field_opt_sf_1A(3,:).*1000;
loops_c=coil_layouts(single_ind_to_plot).out.field_loops_per1Amp(3,:).*1000;
target_c=coil_layouts(single_ind_to_plot).out.target_field_1A.b(3,:).*1000;
target_c = target_c./max(abs(target_c))*max(abs(sf_c)); %scale the target field to the 1A of SF field amplitude for better comparisson
% loops_c_1A=coil_layouts(single_ind_to_plot).out.field_loops_per1Amp(3,:);
% layout_c_1A=coil_layouts(single_ind_to_plot).out.field_layout_per1Amp(3,:);
pos_data=coil_layouts(single_ind_to_plot).out.target_field.coords;

% Calculate the different error metrics
rel_error_sf_target = abs(sf_c-target_c)./max(abs(target_c))*100; %Plot relativ error between target and stream function field
rel_error_layout_target = abs(layout_c-target_c)./max(abs(target_c))*100; %Plot relativ error between layout and target
rel_error_layout_sf = abs(layout_c-sf_c)./max(abs(sf_c))*100; %relativ error between layout and sf field
rel_error_loops_target = abs(loops_c-target_c)./max(abs(target_c))*100; %relativ error between target and unconnected contours
rel_error_loops_layout = abs(loops_c-layout_c)./max(abs(layout_c))*100; %Field difference unconnected contours and final layout


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
plot_colors=rel_error_sf_target;
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
plot_colors=rel_error_layout_target;
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
plot_colors=rel_error_layout_sf;
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
plot_colors=rel_error_loops_target;
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
plot_colors=rel_error_loops_layout;
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

