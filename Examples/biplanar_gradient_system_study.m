
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates a "x" gradient on a biplanar rectangular support
%strcuture

clc; clear all; 

if ispc
cd('..\');
else
cd('../');
end

%% Run the algorithm


%Define the different parameters, which will be varied
sampling_parameters.tikonov_values=10.^(-5:0.1:5);
sampling_parameters.tikonov_values=sampling_parameters.tikonov_values(60);
sampling_parameters.plate_distances=0.1:0.05:0.5;
sampling_parameters.field_shape_functions={'x' 'y' 'z'};
sampling_parameters.plate_alignements={[1 0 0 0] [1 0 0 pi] [0 0 1 pi]};
sampling_parameters.plate_size=0.1:0.05:0.5;

fieldnames=fieldnames(sampling_parameters);
num_fields=numel(fieldnames);
field_numbers_of_elements=arrayfun(@(x) numel(getfield(sampling_parameters,fieldnames{x})),1:num_fields);
field_para_inds=arrayfun(@(x) 1:field_numbers_of_elements(x),1:numel(field_numbers_of_elements),'UniformOutput',false);


%Create a grid for all possible combinations of input parameters 
% (WARNING: DO NOT INSERT TO MANY VALUES;  IT MIGHT CRASH YOUR COMPUTER)
total_number_of_cases=prod(field_numbers_of_elements);
para_ind_grid=cell(1,num_fields);
[para_ind_grid{1:num_fields}] = ndgrid(field_para_inds{:});

result_parameters.sensitivity=zeros(size(para_ind_grid{1}));
result_parameters.mean_error=zeros(size(para_ind_grid{1}));
result_parameters.max_error=zeros(size(para_ind_grid{1}));
result_parameters.wire_length=zeros(size(para_ind_grid{1}));

for case_ind=1:total_number_of_cases

%Select the case parameters
tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(case_ind));
plate_distance=sampling_parameters.plate_distances(para_ind_grid{2}(case_ind));
field_shape_function=sampling_parameters.field_shape_functions{para_ind_grid{3}(case_ind)};
plate_alignement=sampling_parameters.plate_alignements{para_ind_grid{4}(case_ind)};
plate_size=sampling_parameters.plate_size(para_ind_grid{5}(case_ind));

try
    coil_layouts.out=CoilGen(...
    'field_shape_function',field_shape_function,... % definition of the target field
    'coil_mesh_file','create bi-planary mesh', ...
    'biplanar_mesh_parameter_list',[plate_size plate_size 20 20 plate_alignement(1) plate_alignement(2) plate_alignement(3) plate_alignement(4) 0 0 0 plate_distance],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'min_loop_signifcance',3,...
    'target_region_radius',0.01,...  % in meter
    'use_only_target_mesh_verts',false, ...
    'sf_source_file','none', ...
    'levels',14, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',false, ...
    'interconnection_cut_width',0.005, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'skip_postprocessing',true,...
    'skip_inductance_calculation',false,...
    'tikonov_reg_factor',tikonov_value); %Tikonov regularization factor for the SF optimization

result_parameters.max_error(case_ind)=coil_layouts.out.error_vals.max_rel_error_layout_vs_target;
result_parameters.mean_error(case_ind)=coil_layouts.out.error_vals.mean_rel_error_layout_vs_target;
result_parameters.sensitivity(case_ind)=max([coil_layouts.out.layout_gradient.mean_local_gx coil_layouts.out.layout_gradient.mean_local_gy coil_layouts.out.layout_gradient.mean_local_gz]);
%result_parameters.wire_length(case_ind)=sum([coil_layouts.out.coil_parts(:).coil_length]);
result_parameters.wire_length(case_ind)=sum([coil_layouts.out.coil_parts(:).combined_loop_length]);


catch

end


end



%% Select certain parameter cases and plot results

clear inds plot_parameter
inds(1).match=find(para_ind_grid{1}==find(sampling_parameters.tikonov_values==sampling_parameters.tikonov_values));
inds(end+1).match=find(para_ind_grid{2}==find(sampling_parameters.plate_distances==0.2));
inds(end+1).match=find(para_ind_grid{3}==find(strcmp(sampling_parameters.field_shape_functions,'x')));
inds(end+1).match=find(para_ind_grid{4}==find(arrayfun(@(x) isequal(sampling_parameters.plate_alignements{x},[1 0 0 0]),1:numel(sampling_parameters.plate_alignements))));
% inds(end+1).match=find(para_ind_grid{5}==find(sampling_parameters.plate_size==0.2));

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



%% Plot results

close all;


figure;
hold on;
plot(plot_parameter.plate_size,plot_parameter.sensitivity);
grid on;
axis equal;
xlabel('Plate Size [m]', 'Interpreter', 'none');
ylabel('Sensitivity [mT\m\A]', 'Interpreter', 'none');
hold off;




% coil_name='Coil';
% coil_to_plot=coil_layouts(1).out;
% 
% 
% if ispc
% addpath(strcat(pwd,'\','plotting'));
% else
% addpath(strcat(pwd,'/','plotting'));
% end
% %Chose a even leveled solution for plotting
% solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)));
% single_ind_to_plot= find_even_leveled_solution(coil_layouts);
% plot_error_different_solutions(coil_layouts,single_ind_to_plot,coil_name);
% plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
% plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
% plot_coil_parameters(coil_layouts,coil_name);
% plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
% plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
% plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
% rmpath('plotting');


