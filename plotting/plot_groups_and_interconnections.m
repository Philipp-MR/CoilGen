function  plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,plot_title)

if ~coil_layouts.out.input_data.skip_postprocessing

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
%plot the opening cuts within groups
for group_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups)
for cut_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind))
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).cutshape(cut_ind).uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).groups(group_ind).cutshape(cut_ind).uv(2,:),'r','LineWidth',2);
end
end

%plot the opening cuts among groups
for group_ind=1:numel(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups)
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(group_ind).cut1(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(group_ind).cut1(2,:),'color','b','LineWidth',2);
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(group_ind).cut2(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).opening_cuts_among_groups(group_ind).cut2(2,:),'color','b','LineWidth',2);
end
plot(coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.uv(1,:),coil_layouts(single_ind_to_plot).out.coil_parts(part_ind).wire_path.uv(2,:),'k','LineWidth',1);
end
axis equal tight;
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'XColor', 'none','YColor','none');
hold off;
end

set(gcf,'color','w');


end



