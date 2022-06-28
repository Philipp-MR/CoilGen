
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

tikonov_values=10.^(-5:0.5:5);


result_parameters(numel(tikonov_values)).sensitivity=[];
result_parameters(numel(tikonov_values)).mean_error=[];
result_parameters(numel(tikonov_values)).max_error=[];
result_parameters(numel(tikonov_values)).wire_length=[];

for case_ind=1:numel(tikonov_values)
try
    coil_layouts.out=CoilGen(...
    'field_shape_function','x',... % definition of the target field
    'coil_mesh_file','create bi-planary mesh', ...
    'biplanar_mesh_parameter_list',[0.25 0.25 20 20 1 0 0 0 0 0 0 0.25],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
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
    'tikonov_reg_factor',tikonov_values(case_ind)); %Tikonov regularization factor for the SF optimization

result_parameters(case_ind).max_error=coil_layouts.out.error_vals.max_rel_error_layout_vs_target;
result_parameters(case_ind).mean_error=coil_layouts.out.error_vals.mean_rel_error_layout_vs_target;
result_parameters(case_ind).sensitivity=max([coil_layouts.out.layout_gradient.mean_local_gx coil_layouts.out.layout_gradient.mean_local_gy coil_layouts.out.layout_gradient.mean_local_gz]);
%result_parameters(case_ind).wire_length=sum([coil_layouts.out.coil_parts(:).coil_length]);
result_parameters(case_ind).wire_length=sum([coil_layouts.out.coil_parts(:).combined_loop_length]);


catch

end


end


% % % % %% Plot results
% % % % close all;
% % % % 
% % % % coil_name='Coil';
% % % % 
% % % % if ispc
% % % % addpath(strcat(pwd,'\','plotting'));
% % % % else
% % % % addpath(strcat(pwd,'/','plotting'));
% % % % end
% % % % %Chose a even leveled solution for plotting
% % % % solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)));
% % % % single_ind_to_plot= find_even_leveled_solution(coil_layouts);
% % % % plot_error_different_solutions(coil_layouts,single_ind_to_plot,coil_name);
% % % % plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
% % % % plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
% % % % plot_coil_parameters(coil_layouts,coil_name);
% % % % plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
% % % % plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
% % % % plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
% % % % rmpath('plotting');

