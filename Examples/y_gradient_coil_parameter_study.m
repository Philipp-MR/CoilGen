
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates multiple cylindrical "y" gradient coils while
%varying the number of potentials "n" (~the number of turns) 
%All resulting coil parameters for the different solutions are presented in
%the plots


clc; clear all; 

cd('..\');


level_potentials=8:14;

%% Run the algorithm

coil_layouts(numel(level_potentials)).out=[]; %Inititialize the output struct

coil_ind=1;
for level=level_potentials % run the algorithm for different numbers of potential levels

%  try
    coil_layouts(coil_ind).out=CoilGen(...
    'field_shape_function','y',...
    'coil_mesh_file','cylinder_radius500mm_length1500mm.stl', ...    
    'target_mesh_file','none', ... 
    'target_region_radius',0.15,...
    'levels',level, ...
    'pot_offset_factor',0.25, ...
    'surface_is_cylinder_flag',true, ...
    'interconnection_cut_width',0.1, ...
    'normal_shift_length',0.02, ...
    'iteration_num_stream_func_refinement',1, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',...
    'skip_inductance_calculation',false,...
    'tikonov_reg_factor',10);
% catch
% end

coil_ind=coil_ind+1;

end



%% Plot results
close all;

coil_name='coil';

addpath(strcat(pwd,'\','plotting'));
%Chose a even leveled solution for plotting
solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)));
single_ind_to_plot= find_even_leveled_solution(coil_layouts);
plot_error_different_solutions(coil_layouts,solutions_to_plot,coil_name);
plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
%plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
plot_coil_parameters(coil_layouts,coil_name);
plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
rmpath('plotting');

