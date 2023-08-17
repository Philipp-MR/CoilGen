
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This genearets a diffusion weighting MR gradient coil for the female breast. An already
%optimized solution for the stream function is loaded

%For the background of this project refer to: Jia F, Littin S, Amrein P, Yu H, Magill AW, Kuder TA, Bickelhaupt S, Laun F, Ladd ME, Zaitsev M. Design of a high-performance non-linear gradient coil for diffusion weighted MRI of the breast. J Magn Reson. 2021 Oct;331:107052. doi: 10.1016/j.jmr.2021.107052. Epub 2021 Aug 14. PMID: 34478997.


clc; clear all; 

if ispc
cd('..\');
else
cd('../');
end



%% Run the algorithm


    coil_layouts.out=CoilGen(...
    'field_shape_function','none',... % definition of the target field
    'coil_mesh_file','none', ...    
    'min_loop_signifcance',4,...
    'use_only_target_mesh_verts',false, ...
    'sf_source_file','source_data_breast_coil.mat', ...
    'levels',14, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.25, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',false, ...
    'interconnection_cut_width',0.01, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',false,...
    'cross_sectional_points',[sin(0:(2*pi)/(10-1):2*pi); cos(0:(2*pi)/(10-1):2*pi)].*[0.002; 0.002],...
    'skip_sweep',false,...
    'skip_inductance_calculation',false); 




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
plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
plot_coil_parameters(coil_layouts,coil_name);
plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
rmpath('plotting');

