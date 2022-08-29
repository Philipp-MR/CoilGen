function  plot_2D_contours_with_sf(coil_layout,single_ind_to_plot,plot_title)



% Plot a single solution with all steps
if isfield(coil_layout(single_ind_to_plot).out.coil_parts,'groups')
figure('name',plot_title);
tiledlayout('flow');
for part_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts)
nexttile;
title(plot_title+": "+"SF Part"+num2str(part_ind), 'interpreter', 'none');
%Plot the ungrouped and unconnected contour lines with the potential value
%as the color map
hold on;
patch('Faces',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces','Vertices',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.uv','FaceVertexCData',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).stream_function,'FaceColor','interp','edgealpha',0.0);
my_color=parula(numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).potential_level_list));
for group_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups)
for loop_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops)
p=plot(   coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(1,:),...
                coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(2,:),...
                '-o','LineWidth',2,'MarkerSize',0.5);
%color_ind=find(coil_layout(single_ind).out.potential_level_list./coil_layout(single_ind).out.grouping(group_ind).loops(loop_ind).potential==1);
% p.MarkerFaceColor=my_color(color_ind,:);
% p.Color=my_color(color_ind,:);
p.MarkerFaceColor=[0;0;0];
p.Color=[0;0;0];
end
end
%triplot(triangulation(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.uv'),'color',[0.8 0.8 0.8]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis equal tight;
caxis([min(coil_layout(single_ind_to_plot).out.combined_mesh.stream_function) max(coil_layout(single_ind_to_plot).out.combined_mesh.stream_function)]);
set(gca,'XColor', 'none','YColor','none');
set(gcf,'color','w');
hold off;
end
elseif isfield(coil_layout(single_ind_to_plot).out.coil_parts(1),'contour_lines')
figure('name',plot_title);
tiledlayout('flow');
for part_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts)
nexttile;
hold on;
title(plot_title+": "+"SF Part"+num2str(part_ind), 'interpreter', 'none');
patch('Faces',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces','Vertices',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.uv','FaceVertexCData',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).stream_function,'FaceColor','interp','edgealpha',0.0);
for loop_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).contour_lines)
p=plot(   coil_layout.out(single_ind_to_plot).coil_parts(part_ind).contour_lines(loop_ind).uv(1,:),...
              coil_layout.out(single_ind_to_plot).coil_parts(part_ind).contour_lines(loop_ind).uv(2,:),...
                'k-o','LineWidth',1,'MarkerSize',0.5);
end
%triplot(triangulation(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.uv'),'color',[0.8 0.8 0.8]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis equal tight;
caxis([min(coil_layout(single_ind_to_plot).out.combined_mesh.stream_function) max(coil_layout(single_ind_to_plot).out.combined_mesh.stream_function)]);
set(gca,'XColor', 'none','YColor','none');
set(gcf,'color','w');
hold off;
end
end


end

