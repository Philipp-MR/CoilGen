
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates a cylindrical "x" gradient coil for a halbach i.e.
%the cylinder rotated for 90deg 


clc; clear all; 

if ispc
cd('..\');
else
cd('../');
end


%close all;


%% Run the algorithm


    coil_layouts.out=CoilGen(...
    'field_shape_function','x',... % definition of the target field
    'coil_mesh_file','cylinder_radius500mm_length1500mm_90deg_rotated.stl', ...    
    'target_mesh_file','none', ... 
    'secondary_target_mesh_file','none', ...
    'secondary_target_weight',0.5, ...
    'target_region_radius',0.10,...  % in meter
    'use_only_target_mesh_verts',false, ...
    'sf_source_file','none', ...
    'levels',16, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',true, ...
    'interconnection_cut_width',0.035, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'iteration_num_stream_func_refinement',1, ... % the number of refinements for the mesh;
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',false,...
    'skip_inductance_calculation',false,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',1); %Tikonov regularization factor for the SF optimization




%% Plot results
close all;

coil_name='Coil';

if ispc
addpath(strcat(pwd,'\','plotting'));
else
addpath(strcat(pwd,'/','plotting'));
end
%Chose a even leveled solution for plotting
%solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)));
single_ind_to_plot= find_even_leveled_solution(coil_layouts);
plot_error_different_solutions(coil_layouts,single_ind_to_plot,coil_name);
plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
plot_3D_sf(coil_layouts,single_ind_to_plot,coil_name);
plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
plot_coil_parameters(coil_layouts,coil_name);
plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
rmpath('plotting');

