function  plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,plot_title)

if exist('coil_layout(single_ind_to_plot).out.coil_parts(part_ind).groups','var')

figure('name',plot_title);
tiledlayout('flow');
%Plot grouped loops together with the interconnection areas
for part_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts)
nexttile;
hold on;
my_color=jet(numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups));
axis equal tight
for group_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups)
for loop_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops)
p=plot(   coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(1,:),...
                coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(2,:),...
            '-o','LineWidth',4,'MarkerFaceColor',my_color(group_ind,:),'MarkerSize',1.5);
p.Color=my_color(group_ind,:);
end
end
%plot the opening cuts between groups
for cut_ind=1: numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).rectangle_cuts.high)
if  numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(cut_ind).loops)~=1
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).rectangle_cuts.high(cut_ind).uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).rectangle_cuts.high(cut_ind).uv(2,:),'r','LineWidth',3);
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).rectangle_cuts.low(cut_ind).uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).rectangle_cuts.low(cut_ind).uv(2,:),'g','LineWidth',3);
end
end
%plot the opening cuts among groups
for cut_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups)
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(cut_ind).uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(cut_ind).uv(2,:),'color','b','LineWidth',3);
end
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.uv(2,:),'k');
axis equal tight;
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'XColor', 'none','YColor','none');
hold off;
end

end


end

