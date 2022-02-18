function result_out=CoilGen(varargin)


%%%% Create optimized coil finished coil layout 

%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics

%%%% 5.10.2021

% the following external functions were used in modifcated form:
% intreparc@John D'Errico (2010), @matlabcentral/fileexchange
%The non-cylidnrical parameterization is taken from "matlabmesh @ Ryan
%Schmidt  rms@dgp.toronto.edu" based on desbrun et al (2002), "Intrinsic Parameterizations of {Surface} Meshes",
%NS (2021). Curve intersections (https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections), MATLAB Central File Exchange. Retrieved October 6, 2021. 


 

%%% ALGORITHM %%%



addpath(strcat(pwd,'\','sub_functions'));

%Parse the input variables
tic;
disp('Parse inputs:');
[input_parser,input] = parse_input(varargin);
toc; 

if strcmp(input.sf_source_file,'none')

%Read the input mesh
tic;
disp('Load geometry:');
[coil_mesh,target_mesh,secondary_target_mesh]=read_mesh(input);
toc; 

%Split the mesh and the stream function into disconnected pieces
tic;
disp('Split the mesh and the stream function into disconnected pieces.');
coil_parts= split_disconnected_mesh(coil_mesh);
toc; 

%upsample the mesh density by subdivision
tic;
disp('Upsample the mesh by subdivision:');
coil_parts=refine_mesh(coil_parts,input.iteration_num_stream_func_refinement,input.sf_source_file);
toc; 

%Parameterize the mesh
tic;
disp('Parameterize the mesh:');
coil_parts=parameterize_mesh(coil_parts,input.surface_is_cylinder_flag); 
toc; 

%Define the target field
tic;
disp('Define the target field:');
[target_field,is_supressed_point] = define_target_field(coil_parts,target_mesh,secondary_target_mesh,input);
toc; 

%find indices of mesh nodes for one ring basis functions
tic;
disp('Calculate mesh one ring:');
coil_parts=calculate_one_ring_by_mesh(coil_parts);
toc; 

% create the basis funtion container which represents the current density
tic;
disp('Create the basis funtion container which represents the current density:');
coil_parts=calculate_basis_functions(coil_parts);
toc; 

% calculate the sensitivity matrix Cn
tic;
disp('Calculate the sensitivity matrix:');
coil_parts=calculate_sensitivity_matrix(coil_parts,target_field.coords,input.gauss_order);
toc; 

% Calculate the resistance matrix Rmn
tic;
disp('Calculate the resistance matrix:');
coil_parts=calculate_resistance_matrix(coil_parts,input.gauss_order,input.conductor_thickness,input.specific_conductivity_conductor);
toc; 

%Optimize the stream function toward target field and further constraints
tic;
disp('Optimize the stream function toward target field and secondary constraints:');
[coil_parts,combined_mesh,sf_b_field]= stream_function_optimization(coil_parts,target_field,input.tikonov_reg_factor,input.sf_source_file);
toc; 

else % load the preoptimized data
tic;
disp('Load preoptimized data:');
[coil_parts,~,~,combined_mesh,sf_b_field,target_field,is_supressed_point]=load_preoptimized_data(input);
toc; 

end


%Calculate the potential levels for the discretization
tic;
disp('Calculate the potential levels for the discretization:');
[coil_parts,primary_surface_ind]=calc_potential_levels(coil_parts,combined_mesh,input.levels,input.pot_offset_factor,input.level_set_method);
toc; 

%Generate the contours
tic;
disp('Generate the contours:');
coil_parts=calc_contours_by_triangular_potential_cuts(coil_parts);
toc; 

%Process contours
tic;
disp('Process contours:');
coil_parts= process_raw_loops(coil_parts,input.min_point_loop_number,input.area_perimeter_deletion_ratio,input.max_allowed_angle_within_coil_track,input.min_allowed_angle_within_coil_track,input.tiny_segment_length_percentage);
toc; 

if ~input.skip_postprocessing 


%Find the minimal distance between the contour lines
tic;
disp('Find the minimal distance between the contour lines:');
coil_parts= find_minimal_contour_distance(coil_parts,input.track_width_factor);
toc; 

%Group the contour loops in topological order
tic;
disp('Group the contour loops in topological order:');
coil_parts=topological_loop_grouping(coil_parts);
toc; 

%Calculate center locations of groups
tic;
disp('Calculate center locations of groups:');
coil_parts=calculate_group_centers(coil_parts);
toc; 

%open the groups along the b0 direction and surface openings
tic;
disp('Find the interconnection areas:');
coil_parts = calculate_opening_directions(coil_parts,input.b_0_direction,input.interconnection_cut_width);
toc; 


%interconnect the single groups
tic;
disp('Interconnect the single groups:');
coil_parts = interconnect_within_groups(coil_parts,input.force_cut_selection);
toc; 

if ~strcmp(input.interconnection_method,'regular')


tic;
disp('Interconnect the differnt groups: (Spiral In/Out)');
coil_parts = perform_spiral_in_out_intergroup_connection(coil_parts,input);
toc;

    
else % do the spiral in-out double layer interconnection method
    

%interconnect the groups to a single wire path
tic;
disp('Interconnect the groups to a single wire path:');
coil_parts=interconnect_among_groups(coil_parts,input.interconnection_cut_width);
toc; 

%connect the groups and shift the return paths over the surface
tic;
disp('Shift the return paths over the surface:');
%coil_parts = shift_return_paths(coil_parts,input.normal_shift_length);
toc; 
    
end

%Calculate the inductance with fast henry
tic;
disp('Calculate the inductance with fast henry:');
coil_parts = calculate_inductance_by_coil_layout(coil_parts,input.conductor_cross_section_width,input.conductor_cross_section_height,input.skip_inductance_calculation); % output in ohm,Henry,m,m²
toc; 

%Evaluate the result for the final wire track
tic;
disp('Evaluate the result for the final wire track:');
[coil_parts,field_by_layout,field_by_unconnected_loops,field_layout_per1Amp,field_loops_per1Amp,field_error_vals,opt_current_layout] = evaluate_field_errors(coil_parts,target_field,sf_b_field);
toc; 

%Calculate the resuting gradient field
tic;
disp('Calculate the resuting gradient field:');
layout_gradient = calculate_gradient(coil_parts,target_field,field_layout_per1Amp,field_loops_per1Amp,combined_mesh);
toc; 


else
    
    
%Evaluate the result only for the contours
tic;
disp('Evaluate the result for the contours:');
[coil_parts,field_by_layout,field_by_unconnected_loops,field_layout_per1Amp,field_loops_per1Amp,field_error_vals,opt_current_layout,layout_gradient] = evaluate_loop_errors(coil_parts,target_field,sf_b_field);
toc; 
coil_parts(numel(coil_parts)).wire_path.uv=[];
coil_parts(numel(coil_parts)).wire_path.v=[];
end

%Assign the outputs
result_out.coil_parts=coil_parts;
result_out.num_levels=input.levels;
result_out.potential_step=coil_parts(primary_surface_ind).contour_step;
result_out.primary_surface=primary_surface_ind;
result_out.target_field=target_field;
result_out.is_supressed_point=is_supressed_point;
result_out.field_by_layout=field_by_layout;
result_out.field_by_unconnected_loops=field_by_unconnected_loops;
result_out.field_layout_per1Amp=field_layout_per1Amp;
result_out.field_loops_per1Amp=field_loops_per1Amp;
result_out.needed_current_layout=opt_current_layout;
result_out.b_field_opt_sf=sf_b_field;
result_out.error_vals=field_error_vals;
result_out.input_data=input_parser.Results;
result_out.layout_gradient=layout_gradient;
result_out.combined_mesh=combined_mesh;

rmpath('sub_functions');
                        
%%% END OF ALGORITHM %%%




end