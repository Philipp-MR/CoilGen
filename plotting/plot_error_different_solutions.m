function  plot_error_different_solutions(coil_layouts,solutions_to_plot,plot_title)


%valid_layouts=arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)); 
%pot_level=arrayfun(@(x) coil_layouts(x).out.input_data.levels,1:numel(coil_layouts));


% Plot the different errors vs level number
figure('name',plot_title);
hold on;
title({"Max relative field Error [%] vs n" ,"Example"+" "+plot_title}, 'interpreter', 'none');



max_rel_error_layout_vs_target=arrayfun(@(x) coil_layouts(x).out.error_vals.max_rel_error_layout_vs_target,solutions_to_plot);
mean_rel_error_layout_vs_target=arrayfun(@(x) coil_layouts(x).out.error_vals.mean_rel_error_layout_vs_target,solutions_to_plot);
max_rel_error_loops_vs_target=arrayfun(@(x) coil_layouts(x).out.error_vals.max_rel_error_unconnected_contours_vs_target,solutions_to_plot);
mean_rel_error_loops_vs_target=arrayfun(@(x) coil_layouts(x).out.error_vals.mean_rel_error_unconnected_contours_vs_target,solutions_to_plot);
max_rel_error_layout_vs_sf=arrayfun(@(x) coil_layouts(x).out.error_vals.max_rel_error_layout_vs_stream_function_field,solutions_to_plot);
mean_rel_error_layout_vs_sf=arrayfun(@(x) coil_layouts(x).out.error_vals.mean_rel_error_layout_vs_stream_function_field,solutions_to_plot);
max_rel_error_loops_vs_sf=arrayfun(@(x) coil_layouts(x).out.error_vals.max_rel_error_unconnected_contours_vs_stream_function_field,solutions_to_plot);
mean_rel_error_loops_vs_sf=arrayfun(@(x) coil_layouts(x).out.error_vals.mean_rel_error_unconnected_contours_vs_stream_function_field,solutions_to_plot);



%Plot also the error compared to the sf field
p1=plot(solutions_to_plot,max_rel_error_loops_vs_sf,'o-b','linewidth',2);
p2=plot(solutions_to_plot,max_rel_error_layout_vs_sf,'o-r','linewidth',2);
p3=plot(solutions_to_plot,mean_rel_error_loops_vs_sf,'*-b','linewidth',2);
p4=plot(solutions_to_plot,mean_rel_error_layout_vs_sf,'*-r','linewidth',2);
p5=plot(solutions_to_plot,max_rel_error_loops_vs_target,'o-y','linewidth',2);
p6=plot(solutions_to_plot,max_rel_error_layout_vs_target,'o-m','linewidth',2);
p7=plot(solutions_to_plot,mean_rel_error_loops_vs_target,'*-','linewidth',2,'color',[0 0.5 0],'Markerfacecolor',[0 0.5 0]);
p8=plot(solutions_to_plot,mean_rel_error_layout_vs_target,'*-c','linewidth',2);
% p3=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),max_rel_error_loops_vs_sf(mod(solutions_to_plot,2)==0),'--bo','linewidth',3,'Markerfacecolor','b');
% p4=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),max_rel_error_layout_vs_sf(mod(solutions_to_plot,2)==0),':rd','linewidth',3,'Markerfacecolor','r');
%Plot also the error compared to the ideal field
% p5=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),max_rel_error_loops_vs_target(mod(solutions_to_plot,2)==0),'-y','linewidth',3);
% p7=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),mean_rel_error_loops_vs_target(mod(solutions_to_plot,2)==0),'-m','linewidth',3);
% p10=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),max_rel_error_layout_vs_target(mod(solutions_to_plot,2)==0),'-s','linewidth',3,'color',[0 0.5 0]);
% p11=plot(solutions_to_plot(mod(solutions_to_plot,2)==0),mean_rel_error_layout_vs_target(mod(solutions_to_plot,2)==0),'-c','linewidth',3);
% p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; %Select the data to include into the legend
% p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p5.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p6.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p8.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p9.Annotation.LegendInformation.IconDisplayStyle = 'off';
%p10.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p11.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p12.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(arrayfun(@(x) coil_layout(x).out.num_levels,valid_layouts),arrayfun(@(x) coil_layout(x).out.relative_mean_field_error_unconnected_loops,valid_layouts),'-xb');
% plot(arrayfun(@(x) coil_layout(x).out.num_levels,valid_layouts),arrayfun(@(x) coil_layout(x).out.relative_mean_field_error,valid_layouts),'-xr');
grid on;
% if isfield(sf_potential,'numeric_target_field')



legend('Max. Err. [%], Loop-field vs SF-field',...
    'Max. Err. [%], Layout-field vs SF-field'...
    ,'Mean. Err. [%], Loop-field vs SF-field'...
    ,'Mean. Err. [%], Layout-field vs SF-field'...
    ,'Max. Err. [%], Loop-field vs ideal target',...
    'Max. Err. [%], Layout-field vs ideal target'...
    ,'Mean. Err. [%], Loop-field vs ideal target'...
    ,'Mean. Err. [%], Layout-field vs ideal target','FontSize',8);
% else
% legend('Err [%], Unconnected loops','Err [%], Connected wire track','FontSize',13);
% end
%ylim([0 3*mean(error_layout(2:end))]);

ylim([0 max([max_rel_error_loops_vs_sf,max_rel_error_layout_vs_sf,mean_rel_error_loops_vs_sf,mean_rel_error_layout_vs_sf,max_rel_error_loops_vs_target,max_rel_error_layout_vs_target,mean_rel_error_loops_vs_target,mean_rel_error_layout_vs_target])*1.1]);

xticks(solutions_to_plot);


%set(gcf, 'Position',  [1000, 100, 500, 500]);
set(gcf,'color','w');
set(gca,'FontWeight','bold');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a','fontsize',15);
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a','fontsize',15);
xlabel('n','FontSize',15);
ylabel('[%]','FontSize',15);
a = gca; 
a.XAxis.FontSize = 15;
a.YAxis.FontSize = 15;

hold off;



end

