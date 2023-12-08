function [input_parser,input] = parse_input(varargin)



%%%%%%input parameter%%%%
%list of valid input variables


input_parser = inputParser; %create parser object
addParameter(input_parser,'temp',[],@isstruct);
%Add the mesh file that represents the boundary of the target geometry
addParameter(input_parser,'coil_mesh_file','none',@ischar);
%Add the spatial function that defines the field
addParameter(input_parser,'field_shape_function','x',@ischar);
%Define a targeted gradient strength in T/m; (will typically not be
%reached..)
addParameter(input_parser,'target_gradient_strength',1,@isnumeric);
%offset factor for contour levels
addParameter(input_parser,'pot_offset_factor',1/2,@isnumeric);
%file of the target surface mesh, which will define the target field
addParameter(input_parser,'target_mesh_file','none',@ischar);
%shielded_geometry_file; where the field should be supressed
addParameter(input_parser,'secondary_target_mesh_file','none',@ischar);
%weight for the secondary target points 
addParameter(input_parser,'secondary_target_weight',1,@isnumeric);
% flag to use only the target mesh vertices as target coordinates
addParameter(input_parser,'use_only_target_mesh_verts',false,@islogical);
%file of an already optimized stream function
addParameter(input_parser,'sf_source_file','none',@ischar);
%target field file; pointwise definition
addParameter(input_parser,'target_field_definition_file','none',@ischar);
%target field file; fieldname definition
addParameter(input_parser,'target_field_definition_field_name','none',@ischar);
%Stream function optimization method; tikkonov, fmincon
addParameter(input_parser,'sf_opt_method','tikkonov',@ischar);
% %Tikonov regularization factor for the SF optimization
addParameter(input_parser,'tikonov_reg_factor',1,@isnumeric);
% %Parameter for the iterative optimization with fmincon: Number
% iterations, number evalutations, Opti.Toleranze
addParameter(input_parser,'fmincon_parameter',[500 10^10 1.000000e-10 1.000000e-10 1.000000e-10],@isnumeric);
%Number of potential levels
addParameter(input_parser,'levels',10,@isnumeric);
%Specify one of the three ways the level sets are calculated:
%"primary","combined", or "independent"
addParameter(input_parser,'level_set_method','primary',@ischar);
%fieldtype to evaluate; 'gradient' or 'field'
addParameter(input_parser,'fieldtype_to_evaluate',"field",@isstring);
% flag for cylindrical surface; in case of cylinder, a special parameterization will be used
addParameter(input_parser,'surface_is_cylinder_flag',true,@islogical);
%for the cylinder parameterization the ration of outer and inner boundary
addParameter(input_parser,'circular_diameter_factor_cylinder_parameterization',1,@isnumeric);
%the width in meter of the opening cut for the interconnection of the loops
addParameter(input_parser,'interconnection_cut_width',0.01,@isnumeric);
%Radius of a spherical target field
addParameter(input_parser,'target_region_radius',0.15,@isnumeric);
%Number of target points per dimension within the target region
addParameter(input_parser,'target_region_resolution',10,@isnumeric);
% the distance in meter for which crossing lines will be seperated along the normal direction of the surface
addParameter(input_parser,'normal_shift_length',0.001,@isnumeric);
% the minimal required number of point of a single loop; otherwise loops
% will be removed..
addParameter(input_parser,'min_point_loop_number',20,@isnumeric);
%The minimal required field contribution (in percent) to the target field;
%loops that contribute less than that can be deleted
addParameter(input_parser,'min_loop_signifcance',0,@isnumeric);
%additional loop removal criteria which relates to the perimeter to surface
%ratio of the loop
addParameter(input_parser,'area_perimeter_deletion_ratio',5,@isnumeric);
%max allowed angle of the track of the contours
addParameter(input_parser,'max_allowed_angle_within_coil_track',120,@isnumeric);
%min allowed angle of the track of the contours; smaller angles will bef
%converted to straight lines in order to reduce the number of points
addParameter(input_parser,'min_allowed_angle_within_coil_track',0.0001,@isnumeric);
%minimum relative percentage for which points will be deleted which contribute to segments which is extremly short
addParameter(input_parser,'tiny_segment_length_percentage',0,@isnumeric);
%number of refinement iterations of the mesh together with the stream function
addParameter(input_parser,'iteration_num_mesh_refinement',0,@isnumeric);
%The direction (vector) along the interconnections will be aligned
addParameter(input_parser,'b_0_direction',[0;0;1],@isnumeric);
%directory of the .stl  geometry files
if ispc
addParameter(input_parser,'geometry_source_path',strcat(pwd,'\','Geometry_Data'),@ischar);
else
addParameter(input_parser,'geometry_source_path',strcat(pwd,'/','Geometry_Data'),@ischar);
end
%output directory
if ispc
%addParameter(input_parser,'output_directory',strcat(pwd,'\','Results'),@ischar);
addParameter(input_parser,'output_directory',pwd,@ischar);
else
%addParameter(input_parser,'output_directory',strcat(pwd,'/','Results'),@ischar);
addParameter(input_parser,'output_directory',pwd,@ischar);
end
%Flag if the track should be smoothed
addParameter(input_parser,'smooth_flag',true,@islogical);
%smoothing parameter
addParameter(input_parser,'smooth_factor',1,@isnumeric);
%flag to save sweeped .stl
addParameter(input_parser,'save_stl_flag',true,@islogical);
%flag to plot results
addParameter(input_parser,'plot_flag',true,@islogical);
%interconnection_method: Regular or spiral in/out
addParameter(input_parser,'interconnection_method','regular',@ischar);
%Group interconnection_method: 'straight' or 'crossed'
addParameter(input_parser,'group_interconnection_method','crossed',@ischar);
%Flag to skip calculation of minimal winding distance
addParameter(input_parser,'skip_calculation_min_winding_distance',true,@islogical);
%Flag to skip post processing
addParameter(input_parser,'skip_postprocessing',false,@islogical);
%Flag to skip inductance_calculation
addParameter(input_parser,'skip_inductance_calculation',false,@islogical);
%Flag to skip the shifting of return paths
addParameter(input_parser,'skip_normal_shift',false,@islogical);
%smoothing parameters regarding the normal shift
addParameter(input_parser,'normal_shift_smooth_factors',[2 3 2],@isnumeric);
%smoothing parameters regarding the normal vector during the surface sweep
addParameter(input_parser,'normal_sweep_smooth_factor',5,@isnumeric);
%Flag to skip the generation of a volumentric (3D) coil body
addParameter(input_parser,'skip_sweep',false,@islogical);
%Flag to generate a rectangular pcb pattern to wrap around a cylinder
addParameter(input_parser,'make_cylndrical_pcb',false,@islogical);
%Flag to generate a rectangular pcb pattern to wrap around a cylinder
addParameter(input_parser,'pcb_interconnection_method','spiral_in_out',@ischar);
%Factor of shifting the open ends of the spirals in order to avoid
%overlapps; in percent
addParameter(input_parser,'pcb_spiral_end_shift_factor',10,@isnumeric);
%force_cut_selection
addParameter(input_parser,'force_cut_selection',{},@iscell);
%Gaus integration order
addParameter(input_parser,'gauss_order',2,@isnumeric);
% flag to set the roi into the geometric center of the mesh
addParameter(input_parser,'set_roi_into_mesh_center',false,@logical);
%In case of pcb layout, specify the track width
addParameter(input_parser,'track_width_factor',0.5,@isnumeric);
%cross_section_width of the conductor (for the inductance calculation) in meter
addParameter(input_parser,'conductor_cross_section_width',0.002,@isnumeric);
%cross_section_width of the conductor (for the inductance calculation) in meter
addParameter(input_parser,'conductor_cross_section_height',0.002,@isnumeric);
%conducter conductiviy
addParameter(input_parser,'specific_conductivity_conductor',0.018*10^(-6),@isnumeric);
% thickness of the sheet current density of within the stream function representation
addParameter(input_parser,'conductor_thickness',0.005,@isnumeric); 
%2D edge points for direct defintion of the cross section of the conductor
%build circular cut shapes
addParameter(input_parser,'cross_sectional_points',[0 0],@isnumeric);
%specify the paramters for the generation of the (default) cylindrical mesh
% => cylinder_height[in m], cylinder_radius[in m], num_circular_divisions,
% num_longitudinal_divisions, rotation_vector_x, rotation_vector_y, rotation_vector_z, rotation_angle [radian]
addParameter(input_parser,'cylinder_mesh_parameter_list',[0.8 0.3 20 20 1 0 0 0],@isnumeric);
%specify the paramters for the generation of the (default) planar mesh
% => planar_height[in m], planar_radius[in m], num_lateral_divisions,
% num_longitudinal_divisions, rotation_vector_x, rotation_vector_y,
% rotation_vector_z, rotation_angle [radian], center_x [m], center_y [m], center_z [m]
%specify the paramters for the generation of a double cone ("diabolo") shaped mesh
% => cylinder_height[in m],inner_cylinder_height [in m], max_cylinder_radius[in m], min_cylinder_radius[in m], num_circular_divisions,
% num_longitudinal_divisions, rotation_vector_x, rotation_vector_y, rotation_vector_z, rotation_angle [radian]
addParameter(input_parser,'double_cone_mesh_parameter_list',[0.8 0.3  0.3 0.1 20 20 1 0 0 0],@isnumeric);
%specify the paramters for the generation of the (default) planar mesh
% => planar_height[in m], planar_radius[in m], num_lateral_divisions,
% num_longitudinal_divisions, rotation_vector_x, rotation_vector_y,
% rotation_vector_z, rotation_angle [radian], center_x [m], center_y [m], center_z [m]
addParameter(input_parser,'planar_mesh_parameter_list',[0.25 0.25 20 20 1 0 0 0 0 0 0],@isnumeric);
%specify the paramters for the generation of the (default) circular mesh
% => radius[in m], radial_divisions, rotation_vector_x, rotation_vector_y, rotation_vector_z, rotation_angle [radian], center_x [m], center_y[m], center_z [m], 
addParameter(input_parser,'circular_mesh_parameter_list',[0.25 20 1 0 0 0 0 0 0],@isnumeric);
%specify the paramters for the generation of the (default) biplanar mesh
% => biplanar_height[in m], biplanar_radius[in m], num_lateral_divisions,
% num_longitudinal_divisions,
% target_normal_x,target_normal_y,target_normal_z, center_x [m], center_y
% [m], center_z [m], plate distance [mm]
addParameter(input_parser,'biplanar_mesh_parameter_list',[0.25 0.25 20 20 1 0 0 0 0 0 0.2],@isnumeric);
%Parse the input arguments
parse(input_parser,varargin{1}{:});
input=input_parser.Results;

end

