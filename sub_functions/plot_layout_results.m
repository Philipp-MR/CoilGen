function plot_layout_results(sf_opt_b_field,target_field,field_by_layout,field_by_unconnected_loops,coil_mesh,contour_lines,wire_path,opt_stream_func,plot_result)


% Plot resulting field and compare with the target field 

if plot_result

comp_factor=1;

dot_size=200;

layout_c=field_by_layout(3,1:comp_factor:end);
sf_c=sf_opt_b_field(3,1:comp_factor:end);
loops_c=field_by_unconnected_loops(3,1:comp_factor:end);
target_c=target_field.b(3,1:comp_factor:end);
pos_data=target_field.coords(:,1:comp_factor:end);

figure; 
tiledlayout('flow');

nexttile;

%Plot target field
hold on;
axis equal;
title('target Bz, [mT/A]');
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
axis equal;
title('Bz calculated by stream function, [mT/A]');
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
axis equal;
title('Layout Bz, [mT/A]');
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
axis equal;
title('Unconnected Contour Bz, [mT/A]');
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
axis equal;
title('Relative SF error, [%]');
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
axis equal;
title('Relative error layout vs. target, [%]');
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

%Plot relativ error between target and unconnected contours
hold on
axis equal;
title('Relative error unconnected contours vs. target, [%]');
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

%Field differenz unconnected contours and final layou
hold on
axis equal;
title('Field differenz unconnected contours and final layout, [%]');
plot_colors=abs(loops_c-layout_c)./max(abs(target_c))*100;
view(45,45);
colorbar;
colormap(parula);
h = colorbar;
ylabel(h, 'Error %');
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off


% % %Plot target field + current carrying surface
% % figure;
% % hold on;
% % axis equal;
% % title('Bz calculated by stream function, [mT/A]');
% % plot_colors=target_field.b(3,:);
% % view(45,45);
% % colormap(parula);
% % h = colorbar;
% % ylabel(h, '[mT/A]');
% % scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
% % caxis([min(target_field.b(3,:)) max(target_field.b(3,:))]);
% % xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
% % trisurf(triangulation(coil_mesh.faces,coil_mesh.vertices'),'facealpha',0.5,'edgecolor','black','facecolor','cyan');
% % hold off



% Plot the streamfunction interpolatet on 2D flat triangular mesh
figure;
hold on;
axis equal;
title('Stream function by optimization');
normed_pot=opt_stream_func;
normed_pot_1=normed_pot(coil_mesh.faces(:,1));
normed_pot_2=normed_pot(coil_mesh.faces(:,2));
normed_pot_3=normed_pot(coil_mesh.faces(:,3));
%color_mat=color_mat./max(color_mat);
u_vals_1=coil_mesh.uv(1,coil_mesh.faces(:,1))';
u_vals_2=coil_mesh.uv(1,coil_mesh.faces(:,2))';
u_vals_3=coil_mesh.uv(1,coil_mesh.faces(:,3))';
v_vals_1=coil_mesh.uv(2,coil_mesh.faces(:,1))';
v_vals_2=coil_mesh.uv(2,coil_mesh.faces(:,2))';
v_vals_3=coil_mesh.uv(2,coil_mesh.faces(:,3))';
fill([u_vals_1 u_vals_2 u_vals_3]',[v_vals_1 v_vals_2 v_vals_3]',[normed_pot_1 normed_pot_2 normed_pot_3]');
scatter(coil_mesh.uv(1,:),coil_mesh.uv(2,:),100*ones(1,size(coil_mesh.uv,2))',opt_stream_func','filled');
for contour_ind=1:numel(contour_lines)
plot(contour_lines(contour_ind).uv(1,:),contour_lines(contour_ind).uv(2,:),'color','k','linewidth',2);
end
%Together with target field
% plot_colors=target_field_single./max(abs(target_field_single));
% scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
hold off;



% Plot the stream function interpolatet on triangular mesh
figure;
hold on;
axis equal;
title('Stream function by optimization');
view(3);
normed_pot=opt_stream_func;
% normed_pot=(opt_stream_func+abs(min(opt_stream_func)))./max(abs((opt_stream_func+min(opt_stream_func))));
% normed_pot=normed_pot./max(normed_pot);
x_vals_1=coil_mesh.vertices(1,coil_mesh.faces(:,1))';
x_vals_2=coil_mesh.vertices(1,coil_mesh.faces(:,2))';
x_vals_3=coil_mesh.vertices(1,coil_mesh.faces(:,3))';
y_vals_1=coil_mesh.vertices(2,coil_mesh.faces(:,1))';
y_vals_2=coil_mesh.vertices(2,coil_mesh.faces(:,2))';
y_vals_3=coil_mesh.vertices(2,coil_mesh.faces(:,3))';
z_vals_1=coil_mesh.vertices(3,coil_mesh.faces(:,1))';
z_vals_2=coil_mesh.vertices(3,coil_mesh.faces(:,2))';
z_vals_3=coil_mesh.vertices(3,coil_mesh.faces(:,3))';
normed_pot_1=normed_pot(coil_mesh.faces(:,1));
normed_pot_2=normed_pot(coil_mesh.faces(:,2));
normed_pot_3=normed_pot(coil_mesh.faces(:,3));
fill3([x_vals_1 x_vals_2 x_vals_3]',[y_vals_1 y_vals_2 y_vals_3]',[z_vals_1 z_vals_2 z_vals_3]',[normed_pot_1 normed_pot_2 normed_pot_3]');
% %Together with target field
% plot_colors=target_field_single./max(abs(target_field_single));
% colormap(parula);
% scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
for contour_ind=1:numel(contour_lines)
plot3(contour_lines(contour_ind).point_coordinates(1,:),contour_lines(contour_ind).point_coordinates(2,:),contour_lines(contour_ind).point_coordinates(3,:),'color','k','linewidth',2);
end
% plot_colors=abs(loops_c-target_c)./max(abs(target_c))*100;
% scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
% caxis([min(plot_colors) max(plot_colors)]);
% colorbar;
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off;


%Plot layout together with error
figure; 
hold on;
title('Layout and relativ field error');
plot_colors=abs(layout_c-target_c)./max(abs(target_c))*100;
caxis([min(plot_colors) max(plot_colors)]);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
colorbar;
plot3(wire_path.v(1,:),wire_path.v(2,:),wire_path.v(3,:),'color','k','linewidth',2);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
axis equal;
view(45,45);
trisurf(triangulation(coil_mesh.faces,coil_mesh.vertices'),'facealpha',0.5,'edgecolor','black','facecolor','cyan');
hold off;



end


end