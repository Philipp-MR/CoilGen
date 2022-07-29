clc; clear all; close all;

%Select result parameters and plot them

if ispc
cd('..\');
else
cd('../');
end


load Results\results_parameter_study.mat





%% Plot results


plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{'all' [0.2] 'x' [1 0 0 0] 'all'},'Plate Size [m]',5);
plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{'all' 'all' 'x' [1 0 0 0] [0.5]},'Plate Distance [m]',2);




function plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,target_parameters,x_label,para_plot_ind)

inds(numel(target_parameters)).match=[];
plot_parameter.max_error=[];
plot_parameter.mean_error=[];
plot_parameter.sensitivity=[];
plot_parameter.wire_length=[];
plot_parameter.tikonov_value=[];
plot_parameter.plate_distances=[];
plot_parameter.plate_size=[];


if ~strcmp(target_parameters{1},'all')
inds(1).match=find(para_ind_grid{1}==find(sampling_parameters.tikonov_values==target_parameters{1}));
else
inds(1).match=find(para_ind_grid{1});
end
if ~strcmp(target_parameters{2},'all')
inds(2).match=find(para_ind_grid{2}==find(sampling_parameters.plate_distances==target_parameters{2}));
else
inds(2).match=find(para_ind_grid{2});
end
if ~strcmp(target_parameters{3},'all')
inds(3).match=find(para_ind_grid{3}==find(strcmp(sampling_parameters.field_shape_functions,target_parameters{3})));
else
inds(3).match=find(para_ind_grid{3});
end
if ~strcmp(target_parameters{4},'all')
inds(4).match=find(para_ind_grid{4}==find(arrayfun(@(x) isequal(sampling_parameters.plate_alignements{x},target_parameters{4}),1:numel(sampling_parameters.plate_alignements))));
else
inds(4).match=find(para_ind_grid{4});
end
if ~strcmp(target_parameters{5},'all')
inds(5).match=find(para_ind_grid{5}==find(sampling_parameters.plate_size==target_parameters{5}));
else
inds(5).match=find(para_ind_grid{5});
end


%find the case indices which match all the criteria
inds_out=inds(1).match;
if numel(inds)>1
for array_ind=2:numel(inds)
        inds_out= intersect(inds_out,inds(array_ind).match);
end
end

%write the paramters for the matching cases
if ~isempty(inds_out)
plot_parameter.max_error=result_parameters.max_error(inds_out);
plot_parameter.mean_error=result_parameters.mean_error(inds_out);
plot_parameter.sensitivity=result_parameters.sensitivity(inds_out);
plot_parameter.wire_length=result_parameters.wire_length(inds_out);
plot_parameter.tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(inds_out));
plot_parameter.plate_distances=sampling_parameters.plate_distances(para_ind_grid{2}(inds_out));
plot_parameter.plate_size=sampling_parameters.plate_size(para_ind_grid{5}(inds_out));
end

switch para_plot_ind
    case 1
x_plot_vals=plot_parameter.tikonov_values;
    case 2
x_plot_vals=plot_parameter.plate_distances;
    case 3
x_plot_vals=plot_parameter.field_shape_functions;
    case 4
x_plot_vals=plot_parameter.plate_alignements;
    case 5
x_plot_vals=plot_parameter.plate_size;
end

figure;
hold on;
plot(x_plot_vals,plot_parameter.sensitivity,'LineWidth',2);
grid on;
%axis equal;
xlabel(x_label, 'Interpreter', 'none');
ylabel('Sensitivity [mT\m\A]', 'Interpreter', 'none');
axis([min(x_plot_vals) max(x_plot_vals) 0 max(plot_parameter.sensitivity)]);
set(gcf,'color','w');
set(gcf, 'Position',  [100, 100, 400, 400])
hold off;


end



