
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

tikonov_factor=10000;
num_levels=30;
pcb_width=0.002;
cut_width=0.025;
normal_shift=0.006;
min_loop_signifcance=3;
%close all;


circular_resolution=10;
conductor_width=0.0015;
cross_sectional_points=[sin(0:(2*pi)/(circular_resolution-1):2*pi); cos(0:(2*pi)/(circular_resolution-1):2*pi)];
cross_sectional_points=cross_sectional_points.*repmat(conductor_width,[2 1]);
normal_shift_smooth_factors=[5 5 5];

%% Run the algorithm
%try
    coil_x.out=CoilGen(...
    'field_shape_function','x',... % definition of the target field
    'coil_mesh_file','create cylinder mesh', ...    
    'cylinder_mesh_parameter_list',[0.4913 0.154 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',min_loop_signifcance,...
    'target_region_radius',0.1,...  % in meter
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',cut_width, ... % the width for the interconnections are interconnected; in meter
    'conductor_cross_section_width',pcb_width,... %width of the generated pcb tracks
    'normal_shift_length',normal_shift, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'make_cylndrical_pcb',true,...
    'skip_inductance_calculation',false,...
    'cross_sectional_points',cross_sectional_points,...
    'normal_shift_smooth_factors',normal_shift_smooth_factors,...
    'smooth_flag',true,...
    'smooth_factor',2,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',tikonov_factor); %Tikonov regularization factor for the SF optimization
% catch
% end

%try
    coil_y.out=CoilGen(...
    'field_shape_function','y',... % definition of the target field
    'coil_mesh_file','create cylinder mesh', ...    
    'cylinder_mesh_parameter_list',[0.3724 0.166 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',min_loop_signifcance,...
    'target_region_radius',0.1,...  % in meter
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',cut_width, ... % the width for the interconnections are interconnected; in meter
    'conductor_cross_section_width',pcb_width,... %width of the generated pcb tracks
    'normal_shift_length',normal_shift, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'make_cylndrical_pcb',true,...
    'skip_inductance_calculation',false,...
    'cross_sectional_points',cross_sectional_points,...
    'normal_shift_smooth_factors',normal_shift_smooth_factors,...
    'smooth_flag',true,...
    'smooth_factor',2,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',tikonov_factor); %Tikonov regularization factor for the SF optimization
% catch
% end

%try
    coil_z.out=CoilGen(...
    'field_shape_function','z',... % definition of the target field
    'coil_mesh_file','create cylinder mesh', ...    
    'cylinder_mesh_parameter_list',[0.37224 0.177 50 50 0 1 0 pi/2],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'surface_is_cylinder_flag',true, ...
    'min_loop_signifcance',min_loop_signifcance,...
    'target_region_radius',0.1,...  % in meter
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'interconnection_cut_width',cut_width, ... % the width for the interconnections are interconnected; in meter
    'conductor_cross_section_width',pcb_width,... %width of the generated pcb tracks
    'normal_shift_length',normal_shift, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'skip_postprocessing',false,...
    'make_cylndrical_pcb',true,...
    'skip_inductance_calculation',false,...
    'cross_sectional_points',cross_sectional_points,...
    'normal_shift_smooth_factors',normal_shift_smooth_factors,...
    'smooth_flag',true,...
    'smooth_factor',2,...
    'save_stl_flag',true,...
    'tikonov_reg_factor',tikonov_factor); %Tikonov regularization factor for the SF optimization
% catch
% end

%% Plot results
close all;

coil_name='Coil';

coils_to_plot=coil_x;

if ispc
addpath(strcat(pwd,'\','plotting'));
else
addpath(strcat(pwd,'/','plotting'));
end
%Chose a even leveled solution for plotting
%solutions_to_plot=find(arrayfun(@(x) ~isempty(coils_to_plot),1:numel(coil_layouts)));
single_ind_to_plot= find_even_leveled_solution(coils_to_plot);
%plot_error_different_solutions(coils_to_plot,single_ind_to_plot,coil_name);
%plot_2D_contours_with_sf(coils_to_plot,single_ind_to_plot,coil_name);
%plot_3D_sf(coils_to_plot,single_ind_to_plot,coil_name);
plot_groups_and_interconnections(coils_to_plot,single_ind_to_plot,coil_name);
%plot_coil_parameters(coils_to_plot,coil_name);
plot_coil_track_with_resulting_bfield(coils_to_plot,single_ind_to_plot,coil_name);
plot_various_error_metrics(coils_to_plot,single_ind_to_plot,coil_name);
plot_resulting_gradient(coils_to_plot,single_ind_to_plot,coil_name);

plot_pcb_layouts(coils_to_plot,single_ind_to_plot,coil_name);

rmpath('plotting');

