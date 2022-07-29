function  plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,plot_title)



layout_c_1A=coil_layouts(single_ind_to_plot).out.field_layout_per1Amp(3,:);


%Plot the final wire track together with the mesh and the layout field
figure('name',plot_title);
tiledlayout('flow');
nexttile;
hold on;
for part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts)
trisurf(triangulation(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces',coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices'),'facecolor','black','facealpha',0.05,'edgecolor','black','edgealpha',0.05);
if isfield(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind),'wire_path')
plot3(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(2,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(3,:),'linewidth',2,'color',[0 0.4470 0.7410]);
else
for loop_ind=1:numel(coil_layouts.out.coil_parts(part_ind).contour_lines) 
plot3(coil_layouts.out.coil_parts(part_ind).contour_lines(loop_ind).v(1,:),coil_layouts.out.coil_parts(part_ind).contour_lines(loop_ind).v(2,:),coil_layouts.out.coil_parts(part_ind).contour_lines(loop_ind).v(3,:),'linewidth',2,'color',[0 0.4470 0.7410]);
end
end
%title("n = "+num2str(coil_layout(single_ind_to_plot).out.num_levels));
end
xlabel("x[m]");
ylabel("y[m]");
zlabel("z[m]");
%target_sf=coil_layout(single_ind_to_plot).out.target_field;
%plot_vals_scatter=coil_layout(single_ind_to_plot).out.field_by_layout(3,:)./mean(abs(coil_layout(single_ind_to_plot).out.field_by_layout(3,:)));
plot_vals_scatter=coil_layouts.out.field_layout_per1Amp(3,:).*1000;
%plot_vals_scatter=coil_layout(single_ind_to_plot).out.field_by_layout(3,:)./max(abs(coil_layout(single_ind_to_plot).out.field_by_layout(3,:)));
%field_by_layout=coil_layout(single_ind_to_plot).out.field_by_layout(3,:);
%plot_vals_scatter=abs((field_by_layout-target_sf))./max(abs(target_sf)).*100;
h=scatter3(coil_layouts(single_ind_to_plot).out.target_field.coords(1,:),coil_layouts(single_ind_to_plot).out.target_field.coords(2,:),coil_layouts(single_ind_to_plot).out.target_field.coords(3,:),[],plot_vals_scatter,'filled');
caxis([min(plot_vals_scatter), max(plot_vals_scatter)]);
my_bar=colorbar;
title( plot_title+": "+'Layout Bz [mT\\A]', 'interpreter', 'none');
%axis equal tight;
axis equal tight;
view(45,45);
set(gcf,'color','w');
set(gca,'FontWeight','bold');
%set(gca,'FontSize',15);
hold off;

%Plot layout together with error
% nexttile;
% hold on;
% title(plot_title+": "+'Layout vs target field, relative field error [%]');
% plot_colors=abs(layout_c-target_c)./max(abs(target_c))*100;
% caxis([min(plot_colors) max(plot_colors)]);
% %scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
% scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),[],plot_colors,'filled');
% colorbar;
% for part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts)
% plot3(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(2,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.v(3,:),'linewidth',2,'color',[0 0.4470 0.7410]);
% trisurf(triangulation(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces',coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices'),'facecolor','black','facealpha',0.05,'edgecolor','black','edgealpha',0.05);
% end
% xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
% view(45,45);
% set(gcf,'color','w');
% %set(gcf, 'Position',  [1000, 100, 1500, 1000]);
% set(gca,'FontWeight','bold');
% set(gca,'FontSize',15);
% axis equal tight;
% hold off;





end

