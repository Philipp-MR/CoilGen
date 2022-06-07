
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
try
    coil_x.out=CoilGen(...
    'field_shape_function','x',... % definition of the target field
    'cylinder_mesh_parameter_list',[0.4913 0.154 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',1,...
    'target_region_radius',0.1,...  % in meter
    'levels',30, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',0.035, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'skip_inductance_calculation',false,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',100000); %Tikonov regularization factor for the SF optimization
catch
end

try
    coil_y.out=CoilGen(...
    'field_shape_function','y',... % definition of the target field
    'cylinder_mesh_parameter_list',[0.3724 0.166 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',1,...
    'target_region_radius',0.1,...  % in meter
    'levels',30, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',0.035, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'skip_inductance_calculation',false,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',100000); %Tikonov regularization factor for the SF optimization
catch
end

try
    coil_z.out=CoilGen(...
    'field_shape_function','z',... % definition of the target field
    'cylinder_mesh_parameter_list',[0.37224 0.177 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',1,...
    'target_region_radius',0.1,...  % in meter
    'levels',30, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',0.035, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'skip_inductance_calculation',false,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',100000); %Tikonov regularization factor for the SF optimization
catch
end

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
single_ind_to_plot= find_even_leveled_solution(coil_x);
%plot_error_different_solutions(coil_layouts,single_ind_to_plot,coil_name);
%plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
plot_3D_sf(coil_x,single_ind_to_plot,coil_name);
plot_groups_and_interconnections(coil_x,single_ind_to_plot,coil_name);
%plot_coil_parameters(coil_layouts,coil_name);
plot_coil_track_with_resulting_bfield(coil_x,single_ind_to_plot,coil_name);
plot_various_error_metrics(coil_x,single_ind_to_plot,coil_name);
plot_resulting_gradient(coil_x,single_ind_to_plot,coil_name);
rmpath('plotting');

