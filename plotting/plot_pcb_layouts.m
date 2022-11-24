function  plot_pcb_layouts(coil_layouts,single_ind_to_plot,plot_title)

if isfield(coil_layouts(single_ind_to_plot).out.coil_parts(1),'pcb_tracks')

figure;
hold on;
if isfield(coil_layouts(single_ind_to_plot).out.coil_parts.pcb_tracks,'upper_layer')
%plot(layout_2d(1,:),layout_2d(2,:),'b');
for part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts)
for group_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.upper_layer.group_layouts)
for wire_part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.upper_layer.group_layouts(group_ind).wire_parts)
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.upper_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track,'FaceColor','r');
end
end
end

if isfield(coil_layouts(single_ind_to_plot).out.coil_parts.pcb_tracks,'lower_layer')
for group_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.lower_layer.group_layouts)
for wire_part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.lower_layer.group_layouts(group_ind).wire_parts)
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).pcb_tracks.lower_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track,'FaceColor','b');
end
end
end

end
% for wire_part_ind=1:numel(wire_part)
% plot(wire_part(wire_part_ind).uv(1,:)+2*pi,wire_part(wire_part_ind).uv(2,:),'r','linewidth',2);
% end
title(plot_title+": "+"PCB Layout", 'interpreter', 'none');
axis equal;
hold off

set(gcf,'color','w');

end

end

