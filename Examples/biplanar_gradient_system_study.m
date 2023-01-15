
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates a "x" gradient on a biplanar rectangular support
%strcuture

clc; clear all; close all;

if ispc
cd('..\');
else
cd('../');
end



%% Define the different parameters, which will be varied


target_region_radius=0.002;
plate_mesh_resolution=40;

sampling_parameters.tikonov_values=10000000;
% sampling_parameters.plate_distances=[0.01 0.02 0.03 0.04 0.05:0.05:0.5];
% sampling_parameters.plate_size=[0.01 0.02 0.03 0.04 0.05:0.05:0.5];
sampling_parameters.plate_distances=[0.01:0.005:0.1];
sampling_parameters.plate_size=[0.01:0.005:0.1];

% sampling_parameters.plate_distances=[0.05:0.01:0.1];
% sampling_parameters.plate_size=[0.05:0.01:0.1];


sampling_parameters.field_shape_functions={'x'};
sampling_parameters.plate_alignements={[0 0 1]};

% sampling_parameters.tikonov_values=10.^(-5:0.1:5);
% sampling_parameters.tikonov_values=sampling_parameters.tikonov_values(60);
% sampling_parameters.tikonov_values=10000;
% sampling_parameters.plate_distances=0.03;
% sampling_parameters.plate_size=0.05;
% sampling_parameters.field_shape_functions={'x' 'y' 'z'};
% sampling_parameters.plate_alignements=cell(1,2);
% sampling_parameters.plate_alignements{1}=[-0.7485 -0.0000 0.6631];
% sampling_parameters.plate_alignements{2}=[0 0 1];



% sampling_parameters.tikonov_values=100;
% sampling_parameters.plate_distances=[0.1];
% sampling_parameters.plate_size=[0.1];
% sampling_parameters.field_shape_functions={'x'};
% sampling_parameters.plate_alignements=[0 0 1];



% %Calculate different orientation for the gradient system
% azimuth_resol=30;
% polar_resol=30;
% z_val=sin((0:(1/(azimuth_resol-1)):1).*(pi/2));
% cut_radia=(1-(z_val).^2).^(1/2);
% lin_points_num=floor(cut_radia.*polar_resol);
% x_vals=[]; y_vals=[]; z_vals=[]; 
% for az_ind=1:azimuth_resol
% lin_incr=0:1/(lin_points_num(az_ind)-1):1;
% x_vals=[x_vals sin(lin_incr.*(2*pi)).*cut_radia(az_ind)];
% y_vals=[y_vals cos(lin_incr.*(2*pi)).*cut_radia(az_ind)];
% z_vals=[z_vals z_val(az_ind).*ones(size(lin_incr))];
% end
% [coordinate_alignment,~,~] = unique([x_vals; y_vals; z_vals]','rows');
% coordinate_alignment=[coordinate_alignment; [0 0 1]]';
% sampling_parameters.field_shape_functions=cell(1,size(coordinate_alignment,2));
% for align_ind=1:size(coordinate_alignment,2)
% sampling_parameters.field_shape_functions{align_ind}= char("("+num2str(coordinate_alignment(1,align_ind))+")"+".*x+"+"("+num2str(coordinate_alignment(2,align_ind))+")"+".*y+"+"("+num2str(coordinate_alignment(3,align_ind))+").*z");
% end



% %% Calculate different orientation for the gradient system
% azimuth_resol=30;
% polar_resol=30;
% z_val=sin((0:(1/(azimuth_resol-1)):1).*(pi/2));
% cut_radia=(1-(z_val).^2).^(1/2);
% lin_points_num=floor(cut_radia.*polar_resol);
% x_vals=[]; y_vals=[]; z_vals=[]; 
% for az_ind=1:azimuth_resol
% lin_incr=0:1/(lin_points_num(az_ind)-1):1;
% x_vals=[x_vals sin(lin_incr.*(2*pi)).*cut_radia(az_ind)];
% y_vals=[y_vals cos(lin_incr.*(2*pi)).*cut_radia(az_ind)];
% z_vals=[z_vals z_val(az_ind).*ones(size(lin_incr))];
% end
% [bi_planar_alignment,~,~] = unique([x_vals; y_vals; z_vals]','rows');
% bi_planar_alignment=[[0 0 1]; bi_planar_alignment]';
% sampling_parameters.plate_alignements=cell(1,size(bi_planar_alignment,2));
% for align_ind=1:size(bi_planar_alignment,2)
% sampling_parameters.plate_alignements{align_ind}=[bi_planar_alignment(:,align_ind)'];
% end


% line_resolution=50;
% line_vals=[(-1):2/line_resolution:0 0:2/line_resolution:1];
% x_vals=[zeros(1,numel(line_vals)) line_vals];
% y_vals=[line_vals zeros(1,numel(line_vals))];
% r_vals=x_vals.^2+y_vals.^2;
% z_vals=(1-r_vals.^2).^(1/2);
% [bi_planar_alignment,~,~] = unique([x_vals; y_vals; z_vals]','rows');
% bi_planar_alignment=[[0 0 1]; bi_planar_alignment]';
% sampling_parameters.plate_alignements=cell(1,size(bi_planar_alignment,2));
% for align_ind=1:size(bi_planar_alignment,2)
% sampling_parameters.plate_alignements{align_ind}=[bi_planar_alignment(:,align_ind)'];
% end


%% Build parameter sets of all possibel parameter combinations

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

%% Run the algorithm for the different parameter sets

%coil_layouts(total_number_of_cases).out=[];
coil_layouts(total_number_of_cases).out=[];

for case_ind=1:total_number_of_cases
%for case_ind=1

%Select the case parameters
tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(case_ind));
plate_distance=sampling_parameters.plate_distances(para_ind_grid{2}(case_ind));
plate_size=sampling_parameters.plate_size(para_ind_grid{3}(case_ind));
field_shape_function=sampling_parameters.field_shape_functions{para_ind_grid{4}(case_ind)};
plate_alignement=sampling_parameters.plate_alignements{para_ind_grid{5}(case_ind)};

%field_shape_function= char("("+num2str(plate_alignement(1))+")"+".*x+"+"("+num2str(plate_alignement(2))+")"+".*y+"+"("+num2str(plate_alignement(3))+").*z");

  %try
    coil_layouts(case_ind).out=CoilGen(...
    'field_shape_function',field_shape_function,... % definition of the target field
    'coil_mesh_file','create bi-planary mesh', ...
    'biplanar_mesh_parameter_list',[plate_size plate_size plate_mesh_resolution plate_mesh_resolution plate_alignement(1) plate_alignement(2) plate_alignement(3) 0 0 0 plate_distance],... % cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,  num_longitudinal_divisions, rotation_vector: x,y,z, and  rotation_angle [radian]
    'min_loop_signifcance',1,...
    'target_region_radius',target_region_radius,...  % in meter
    'use_only_target_mesh_verts',false, ...
    'sf_source_file','none', ...
    'levels',14, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.5, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',false, ...
    'interconnection_cut_width',0.005, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.01, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'level_set_method','combined',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'skip_postprocessing',true,...
    'skip_inductance_calculation',false,...
    'skip_normal_shift',true,...
    'skip_sweep',true,...
    'make_cylndrical_pcb',false,...
    'tikonov_reg_factor',tikonov_value); %Tikonov regularization factor for the SF optimization

% result_parameters.max_error(case_ind)=coil_layouts(case_ind).out.error_vals.max_rel_error_layout_vs_target;
% result_parameters.mean_error(case_ind)=coil_layouts(case_ind).out.error_vals.mean_rel_error_layout_vs_target;

switch field_shape_function
    case 'x'
result_parameters.sensitivity(case_ind)=mean(coil_layouts(case_ind).out.layout_gradient.dBzdxyz(1,:));
    case 'y'
result_parameters.sensitivity(case_ind)=mean(coil_layouts(case_ind).out.layout_gradient.dBzdxyz(2,:));
    case 'z'
result_parameters.sensitivity(case_ind)=mean(coil_layouts(case_ind).out.layout_gradient.dBzdxyz(3,:));
end

%[result_parameters.sensitivity(case_ind),~]=calc_gradient_along_vector(coil_layouts(case_ind).out.field_by_layout,coil_layouts(case_ind).out.target_field.coords,coil_layouts(case_ind).out.input_data.field_shape_function);
result_parameters.wire_length(case_ind)=sum([coil_layouts(case_ind).out.coil_parts(:).combined_loop_length]);
result_parameters.max_error(case_ind)=coil_layouts(case_ind).out.error_vals.max_rel_error_layout_vs_target;
result_parameters.mean_error(case_ind)=coil_layouts(case_ind).out.error_vals.mean_rel_error_layout_vs_target;


if case_ind>1
coil_layouts(case_ind).out=nan;
end

% catch
% 
% end


end



% %% Select certain parameter cases and plot results
% clear inds plot_parameter
% inds(1).match=find(para_ind_grid{1}==find(sampling_parameters.tikonov_values==sampling_parameters.tikonov_values));
% inds(end+1).match=find(para_ind_grid{2}==find(sampling_parameters.plate_distances==0.2));
% inds(end+1).match=find(para_ind_grid{3}==find(strcmp(sampling_parameters.field_shape_functions,'x')));
% %inds(end+1).match=find(para_ind_grid{4}==find(arrayfun(@(x) isequal(sampling_parameters.plate_alignements{x},[1 0 0 0]),1:numel(sampling_parameters.plate_alignements))));
% inds(end+1).match=find(para_ind_grid{5}==find(sampling_parameters.plate_size==0.2));
% 
% %find the case indices which match all the criteria
% inds_out=inds(1).match;
% if numel(inds)>1
% for array_ind=2:numel(inds)
%         inds_out= intersect(inds_out,inds(array_ind).match);
% end
% end
% 
% %write the paramters for the matching cases
% if ~isempty(inds_out)
% plot_parameter.max_error=result_parameters.max_error(inds_out);
% plot_parameter.mean_error=result_parameters.mean_error(inds_out);
% plot_parameter.sensitivity=result_parameters.sensitivity(inds_out);
% plot_parameter.wire_length=result_parameters.wire_length(inds_out);
% plot_parameter.tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(inds_out));
% plot_parameter.plate_distances=sampling_parameters.plate_distances(para_ind_grid{2}(inds_out));
% plot_parameter.plate_size=sampling_parameters.plate_size(para_ind_grid{5}(inds_out));
% end
% 
% 
% figure;
% hold on;
% plot(plot_parameter.plate_size,plot_parameter.sensitivity);
% grid on;
% axis equal;
% xlabel('Plate Size [m]', 'Interpreter', 'none');
% ylabel('Sensitivity [mT\m\A]', 'Interpreter', 'none');
% hold off;


%% Plot results

close all;

coil_name='Coil';
coil_to_plot=coil_layouts(1);


if ispc
addpath(strcat(pwd,'\','plotting'));
else
addpath(strcat(pwd,'/','plotting'));
end
%Chose a even leveled solution for plotting
solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_to_plot),1:numel(coil_to_plot)));
%plot_error_different_solutions(coil_to_plot,1,coil_name);
plot_2D_contours_with_sf(coil_to_plot,1,coil_name);
plot_groups_and_interconnections(coil_to_plot,1,coil_name);
%plot_coil_parameters(coil_to_plot,coil_name);
plot_coil_track_with_resulting_bfield(coil_to_plot,1,coil_name);
plot_various_error_metrics(coil_to_plot,1,coil_name);
plot_resulting_gradient(coil_to_plot,1,coil_name);
rmpath('plotting');





