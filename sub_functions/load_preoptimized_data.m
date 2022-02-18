
function [coil_parts,coil_mesh,secondary_target_mesh,combined_mesh,sf_b_field,target_field,is_supressed_point]=load_preoptimized_data(input)

%load preoptimized data
load(strcat('Pre_Optimized_Solutions\',input.sf_source_file))

%Split the mesh and the stream function into disconnected pieces
tic;
disp('Split the mesh and the stream function into disconnected pieces.');
coil_parts= split_disconnected_mesh(coil_mesh);
toc; 

%Parameterize the mesh
tic;
disp('Parameterize the mesh:');
coil_parts=parameterize_mesh(coil_parts,input.surface_is_cylinder_flag); 
toc; 

secondary_target_mesh=[];
combined_mesh=coil_parts.coil_mesh;
combined_mesh.bounding_box=[min(combined_mesh.vertices(1,:)) max(combined_mesh.vertices(1,:)); min(combined_mesh.vertices(2,:)) max(combined_mesh.vertices(2,:)); min(combined_mesh.vertices(3,:)) max(combined_mesh.vertices(3,:))];
combined_mesh.stream_function=stream_function;


coil_mesh=coil_parts.coil_mesh;
coil_parts.stream_function=stream_function;
target_field.weigths=ones(1,size(target_field.b,2));
target_field.target_field_group_inds=ones(1,size(target_field.b,2));
is_supressed_point=zeros(1,size(target_field.b,2));
sf_b_field=target_field.b;

end
