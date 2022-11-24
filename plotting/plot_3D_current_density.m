function  plot_3D_current_density(coil_layout,single_ind_to_plot,plot_title)

% Plot the stream function interpolatet on triangular mesh

figure;
hold on;
title(plot_title+": "+'Current Density by optimization and target Bz', 'interpreter', 'none');
view(3);
colormap hot;


for part_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts)
	x_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    x_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    x_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    y_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    y_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    y_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    z_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
    z_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
    z_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));

    normed_pot=vecnorm(coil_layout.out.coil_parts(part_ind).current_density)';
    normed_pot=normed_pot./max(abs(normed_pot));
%     normed_pot_1=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
%     normed_pot_2=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
%     normed_pot_3=normed_pot(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
    fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',normed_pot,'EdgeAlpha',0);

%     %Plot also the mesh
%     trisurf(triangulation(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces',coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices'),...
%     'facecolor','black','facealpha',0,'edgecolor','black','edgealpha',0.5,'linewidth',1);

%     %Highlight the mesh vertices
%     num_verts=size(coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices,2);
%     scatter3(...
%                 coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,:),...
%                 coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,:),...
%                 coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,:),...
%                 ones(1,num_verts)*10,repmat([0 0 0],[num_verts 1]),'filled');

    %Plot the Iso-Contour lines
    for loop_ind=1:numel(coil_layout.out.coil_parts(part_ind).contour_lines) 
    plot3(coil_layout.out.coil_parts(part_ind).contour_lines(loop_ind).v(1,:).*1.005,coil_layout.out.coil_parts(part_ind).contour_lines(loop_ind).v(2,:).*1.005,coil_layout.out.coil_parts(part_ind).contour_lines(loop_ind).v(3,:).*1.005,'linewidth',2,'color',[0 0.4470 0.7410]);
    end


end

box off;
plot_out = gca;
plot_out.XAxis.Visible = 'off';
plot_out.YAxis.Visible = 'off';
plot_out.ZAxis.Visible = 'off';
set(gcf,'color','w');
axis equal;
hold off;

end



% % for part_ind=1:numel(coil_layout(single_ind_to_plot).out.coil_parts)
% %     normed_pot=coil_layout(single_ind_to_plot).out.coil_parts.vertex_current_density;
% %     %normed_pot=normed_pot./max(abs(normed_pot));
% % 	x_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
% %     x_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
% %     x_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(1,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
% %     y_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
% %     y_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
% %     y_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(2,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
% %     z_vals_1=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:));
% %     z_vals_2=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:));
% %     z_vals_3=coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.vertices(3,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:));
% %     normed_pot_1=vecnorm(normed_pot(:,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(1,:)));
% %     normed_pot_2=vecnorm(normed_pot(:,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(2,:)));
% %     normed_pot_3=vecnorm(normed_pot(:,coil_layout(single_ind_to_plot).out.coil_parts(part_ind).coil_mesh.faces(3,:)));
% %     fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',[normed_pot_1' normed_pot_2' normed_pot_3']','EdgeAlpha',0);
% % 
% % end

% % % function plot_sf_on_mesh(coil_parts,target_field,optimized_field)
% % % 
% % % % Plot the stream function interpolatet on triangular mesh
% % % figure;
% % % hold on;
% % % h(numel(coil_parts)+1).plot=[];
% % % for part_ind=1:numel(coil_parts)
% % %      
% % % normed_pot=coil_parts(part_ind).stream_function;
% % % normed_pot=(normed_pot-min(normed_pot))./max(abs(normed_pot-min(normed_pot)));
% % % x_vals_1=coil_parts(part_ind).coil_mesh.vertices(1,coil_parts(part_ind).coil_mesh.faces(1,:))';
% % % x_vals_2=coil_parts(part_ind).coil_mesh.vertices(1,coil_parts(part_ind).coil_mesh.faces(2,:))';
% % % x_vals_3=coil_parts(part_ind).coil_mesh.vertices(1,coil_parts(part_ind).coil_mesh.faces(3,:))';
% % % y_vals_1=coil_parts(part_ind).coil_mesh.vertices(2,coil_parts(part_ind).coil_mesh.faces(1,:))';
% % % y_vals_2=coil_parts(part_ind).coil_mesh.vertices(2,coil_parts(part_ind).coil_mesh.faces(2,:))';
% % % y_vals_3=coil_parts(part_ind).coil_mesh.vertices(2,coil_parts(part_ind).coil_mesh.faces(3,:))';
% % % z_vals_1=coil_parts(part_ind).coil_mesh.vertices(3,coil_parts(part_ind).coil_mesh.faces(1,:))';
% % % z_vals_2=coil_parts(part_ind).coil_mesh.vertices(3,coil_parts(part_ind).coil_mesh.faces(2,:))';
% % % z_vals_3=coil_parts(part_ind).coil_mesh.vertices(3,coil_parts(part_ind).coil_mesh.faces(3,:))';
% % % normed_pot_1=normed_pot(coil_parts(part_ind).coil_mesh.faces(1,:));
% % % normed_pot_2=normed_pot(coil_parts(part_ind).coil_mesh.faces(2,:));
% % % normed_pot_3=normed_pot(coil_parts(part_ind).coil_mesh.faces(3,:));
% % % h(part_ind).plot=fill3([x_vals_1 x_vals_2 x_vals_3]',[y_vals_1 y_vals_2 y_vals_3]',[z_vals_1 z_vals_2 z_vals_3]',[normed_pot_1 normed_pot_2 normed_pot_3]');
% % % 
% % % %set(h(part_ind).plot,'caxis',[min(normed_pot) max(normed_pot)]);
% % % 
% % % end
% % % 
% % % 
% % % % %Together with target field
% % % title('Stream function by optimization');
% % % view(3);
% % % target_c=optimized_field(3,:);
% % % target_c=(target_c-min(target_c))./max(abs(target_c-min(target_c)));
% % % pos_data=target_field.coords;
% % % h(numel(coil_parts)+1).plot=scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),50*ones(1,numel(pos_data,2)),target_c,'filled');
% % % %set(h(numel(coil_parts)+1).plot,'caxis',[min(target_c) max(target_c)]);
% % % axis equal tight;
% % % xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
% % % hold off;
% % % 
% % % 
% % % 
% % % end
