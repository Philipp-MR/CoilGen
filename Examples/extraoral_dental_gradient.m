
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates a cylindrical "y" gradient coil


clc; clear all; close all;



if ispc
cd('..\');
else
cd('../');
end


close all;


%% Run the algorithm

  %  'field_shape_function','(x.*x+0.7.*y.*y).^(1/2)',... % definition of the target field

%tik_factor=4000;
%tik_factor=[1000 1500 2000 2500 3000];
tik_factor=[1:10:100 150:50:500 600:100:1000 2000:1000:10000];
%tik_factor=[1:10:100 150:50:500 600:100:1000 2000:1000:10000 20000:10000:50000];
num_iteration=5000;
fmincon_parameter=[num_iteration 10^10 1.000000e-14 10^(-18) 10^(-18)];
num_levels=40;
optim_method='tikkonov';
target_region_res=10;
skip_postprocessing=true;
iteration_num_mesh_refinement=1;
min_loop_signifcance=0.001;
target_field_definition_file='intraoral_dental_target_field.mat';
coil_mesh_file='dental_extraoral_ccs2_shifted_2cm_in_y.stl'; 
surface_is_cylinder_flag=true;
normal_shift_length=0.001;
pot_offset_factor=0.25;
interconnection_cut_width=0.002;
skip_inductance_calculation=false;
temp.init=[]; %this struct must be set for the first time


coil_x(numel(tik_factor)).out=[];
coil_y(numel(tik_factor)).out=[];
coil_z(numel(tik_factor)).out=[];
coil_r(numel(tik_factor)).out=[];
coil_phi(numel(tik_factor)).out=[];

for coil_ind=1:numel(coil_x)
% try
    [coil_x(coil_ind).out,temp]=CoilGen(...
    'temp',temp,...
    'target_field_definition_file',target_field_definition_file,...
    'target_field_definition_field_name','lin_x',...
    'coil_mesh_file',coil_mesh_file, ... 
    'min_loop_signifcance',min_loop_signifcance,...
    'iteration_num_mesh_refinement',iteration_num_mesh_refinement,...
    'target_region_resolution',target_region_res,...
    'use_only_target_mesh_verts',false, ...
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',pot_offset_factor, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',surface_is_cylinder_flag, ...
    'interconnection_cut_width',interconnection_cut_width, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',normal_shift_length, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',skip_postprocessing,...
    'skip_sweep',true,...
    'skip_inductance_calculation',skip_inductance_calculation,...
    'sf_opt_method',optim_method,...
    'fmincon_parameter',fmincon_parameter,...
    'tikonov_reg_factor',tik_factor(coil_ind)); %Tikonov regularization factor for the SF optimization
% catch
% end
coil_x(coil_ind).out.coil_parts.sensitivity_matrix=[];
coil_x(coil_ind).out.coil_parts.gradient_sensitivity_matrix=[];
coil_x(coil_ind).out.coil_parts.current_density_mat=[];
coil_x(coil_ind).out.coil_parts.triangle_corner_coord_mat=[];
coil_x(coil_ind).out.coil_parts.node_triangle_mat=[];
coil_x(coil_ind).out.coil_parts.basis_elements=[];
end


for coil_ind=1:numel(coil_y)
try
    coil_y(coil_ind).out=CoilGen(...
    'target_field_definition_file',target_field_definition_file,...
    'target_field_definition_field_name','lin_y',...
    'coil_mesh_file',coil_mesh_file, ... 
    'min_loop_signifcance',min_loop_signifcance,...
    'iteration_num_mesh_refinement',iteration_num_mesh_refinement,...
    'target_region_resolution',target_region_res,...
    'use_only_target_mesh_verts',false, ...
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',pot_offset_factor, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',surface_is_cylinder_flag, ...
    'interconnection_cut_width',interconnection_cut_width, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',normal_shift_length, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',skip_postprocessing,...
    'skip_sweep',true,...
    'skip_inductance_calculation',skip_inductance_calculation,...
    'sf_opt_method',optim_method,...
    'fmincon_parameter',fmincon_parameter,...
    'tikonov_reg_factor',tik_factor(coil_ind)); %Tikonov regularization factor for the SF optimization
catch
end
coil_y(coil_ind).out.coil_parts.sensitivity_matrix=[];
coil_y(coil_ind).out.coil_parts.gradient_sensitivity_matrix=[];
coil_y(coil_ind).out.coil_parts.current_density_mat=[];
coil_y(coil_ind).out.coil_parts.triangle_corner_coord_mat=[];
coil_y(coil_ind).out.coil_parts.node_triangle_mat=[];
coil_y(coil_ind).out.coil_parts.basis_elements=[];
end


for coil_ind=1:numel(coil_z)
try
    coil_z(coil_ind).out=CoilGen(...
    'target_field_definition_file',target_field_definition_file,...
    'target_field_definition_field_name','lin_z',...
    'coil_mesh_file',coil_mesh_file, ... 
    'min_loop_signifcance',min_loop_signifcance,...
    'iteration_num_mesh_refinement',iteration_num_mesh_refinement,...
    'target_region_resolution',target_region_res,...
    'use_only_target_mesh_verts',false, ...
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',pot_offset_factor, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',surface_is_cylinder_flag, ...
    'interconnection_cut_width',interconnection_cut_width, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',normal_shift_length, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',skip_postprocessing,...
    'skip_sweep',true,...
    'skip_inductance_calculation',skip_inductance_calculation,...
    'sf_opt_method',optim_method,...
    'fmincon_parameter',fmincon_parameter,...
    'tikonov_reg_factor',tik_factor(coil_ind)); %Tikonov regularization factor for the SF optimization
catch
end
coil_z(coil_ind).out.coil_parts.sensitivity_matrix=[];
coil_z(coil_ind).out.coil_parts.gradient_sensitivity_matrix=[];
coil_z(coil_ind).out.coil_parts.current_density_mat=[];
coil_z(coil_ind).out.coil_parts.triangle_corner_coord_mat=[];
coil_z(coil_ind).out.coil_parts.node_triangle_mat=[];
coil_z(coil_ind).out.coil_parts.basis_elements=[];
end


for coil_ind=1:numel(coil_r)
% try
    coil_r(coil_ind).out=CoilGen(...
    'target_field_definition_file',target_field_definition_file,...
    'target_field_definition_field_name','r',...
    'coil_mesh_file',coil_mesh_file, ... 
    'min_loop_signifcance',min_loop_signifcance,...
    'iteration_num_mesh_refinement',iteration_num_mesh_refinement,...
    'target_region_resolution',target_region_res,...
    'use_only_target_mesh_verts',false, ...
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',pot_offset_factor, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',surface_is_cylinder_flag, ...
    'interconnection_cut_width',interconnection_cut_width, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',normal_shift_length, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',skip_postprocessing,...
    'skip_sweep',true,...
    'skip_inductance_calculation',skip_inductance_calculation,...
    'sf_opt_method',optim_method,...
    'fmincon_parameter',fmincon_parameter,...
    'tikonov_reg_factor',tik_factor(coil_ind)); %Tikonov regularization factor for the SF optimization
% catch
% end
coil_r(coil_ind).out.coil_parts.sensitivity_matrix=[];
coil_r(coil_ind).out.coil_parts.gradient_sensitivity_matrix=[];
coil_r(coil_ind).out.coil_parts.current_density_mat=[];
coil_r(coil_ind).out.coil_parts.triangle_corner_coord_mat=[];
coil_r(coil_ind).out.coil_parts.node_triangle_mat=[];
coil_r(coil_ind).out.coil_parts.basis_elements=[];
end


for coil_ind=1:numel(coil_phi)
try
    coil_phi(coil_ind).out=CoilGen(...
    'target_field_definition_file',target_field_definition_file,...
    'target_field_definition_field_name','r',...
    'coil_mesh_file',coil_mesh_file, ... 
    'min_loop_signifcance',min_loop_signifcance,...
    'iteration_num_mesh_refinement',iteration_num_mesh_refinement,...
    'target_region_resolution',target_region_res,...
    'use_only_target_mesh_verts',false, ...
    'levels',num_levels, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',pot_offset_factor, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',surface_is_cylinder_flag, ...
    'interconnection_cut_width',interconnection_cut_width, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',normal_shift_length, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'set_roi_into_mesh_center',true, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',skip_postprocessing,...
    'skip_sweep',true,...
    'skip_inductance_calculation',skip_inductance_calculation,...
    'sf_opt_method',optim_method,...
    'fmincon_parameter',fmincon_parameter,...
    'tikonov_reg_factor',tik_factor(coil_ind)); %Tikonov regularization factor for the SF optimization
catch
end
coil_phi(coil_ind).out.coil_parts.sensitivity_matrix=[];
coil_phi(coil_ind).out.coil_parts.gradient_sensitivity_matrix=[];
coil_phi(coil_ind).out.coil_parts.current_density_mat=[];
coil_phi(coil_ind).out.coil_parts.triangle_corner_coord_mat=[];
coil_phi(coil_ind).out.coil_parts.node_triangle_mat=[];
coil_phi(coil_ind).out.coil_parts.basis_elements=[];
end


%% Plot results
close all;


coil_to_plot=coil_r;


tik_vals=arrayfun(@(x) coil_to_plot(x).out.input_data.tikonov_reg_factor,1:numel(coil_to_plot));
error_vals=arrayfun(@(x) coil_to_plot(x).out.error_vals.max_rel_error_layout_vs_target,1:numel(coil_to_plot));
sensitvity_vals=arrayfun(@(x) coil_to_plot(x).out.layout_gradient.mean_gradient_in_target_direction,1:numel(coil_to_plot));

figure;
hold on;
title('Sensitivity');
plot(tik_vals,sensitvity_vals);
axis([0 max(tik_vals) 0 max(sensitvity_vals)])
grid on;
hold off;


figure;
hold on;
title('error');
plot(tik_vals,error_vals);
axis([0 max(tik_vals) 0 max(error_vals)])
grid on;
hold off;


% figure;
% hold on;
% title('Error^2/Senstivity');
% plot(tik_vals,error_vals.*error_vals./sensitvity_vals);
% axis([0 max(tik_vals) 0 max(error_vals.*error_vals./sensitvity_vals)])
% grid on;
% hold off;


[~,min_ind]=min(error_vals);
coil_to_plot=coil_to_plot(min_ind);
%coil_to_plot=coil_to_plot(1);

if ispc
addpath(strcat(pwd,'\','plotting'));
else
addpath(strcat(pwd,'/','plotting'));
end
coil_name='Coil';
%Chose a even leveled solution for plotting
%solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_to_plot(x).out),1:numel(coil_layouts)));
single_ind_to_plot= find_even_leveled_solution(coil_to_plot);
%plot_error_different_solutions(coil_to_plot,single_ind_to_plot,coil_name);
% plot_2D_contours_with_sf(coil_to_plot,single_ind_to_plot,coil_name);
% plot_3D_sf(coil_to_plot,single_ind_to_plot,coil_name);
plot_groups_and_interconnections(coil_to_plot,single_ind_to_plot,coil_name);
plot_coil_parameters(coil_to_plot,coil_name);
plot_coil_track_with_resulting_bfield(coil_to_plot,single_ind_to_plot,coil_name);
plot_various_error_metrics(coil_to_plot,single_ind_to_plot,coil_name);
std_gradient=plot_resulting_gradient(coil_to_plot,single_ind_to_plot,coil_name);
rmpath('plotting');

