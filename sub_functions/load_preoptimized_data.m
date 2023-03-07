function [coil_parts,coil_mesh,secondary_target_mesh,combined_mesh,sf_b_field,target_field,is_supressed_point]=load_preoptimized_data(input)

%load preoptimized data


if ispc
load(strcat('Pre_Optimized_Solutions\',input.sf_source_file));
else
load(strcat('Pre_Optimized_Solutions/',input.sf_source_file));
end

secondary_target_mesh=[];

%Split the mesh and the stream function into disconnected pieces
tic;
disp('Split the mesh and the stream function into disconnected pieces.');
coil_parts= split_disconnected_mesh(coil_mesh);
toc; 

%Parameterize the mesh
tic;
disp('Parameterize the mesh:');
coil_parts=parameterize_mesh(coil_parts,input); 
toc; 


%define target field
target_field.weigths=ones(1,size(target_field.b,2));
target_field.target_field_group_inds=ones(1,size(target_field.b,2));
is_supressed_point=zeros(1,size(target_field.b,2));
sf_b_field=target_field.b;


%generate a combined mesh container
combined_mesh=coil_mesh;
combined_mesh.bounding_box=[min(combined_mesh.vertices(1,:)) max(combined_mesh.vertices(1,:)); min(combined_mesh.vertices(2,:)) max(combined_mesh.vertices(2,:)); min(combined_mesh.vertices(3,:)) max(combined_mesh.vertices(3,:))];
combined_mesh.stream_function=stream_function;
combined_mesh.faces=coil_parts(1).coil_mesh.faces;
combined_mesh.vertices=coil_parts(1).coil_mesh.vertices;
combined_mesh.n=coil_parts(1).coil_mesh.n;
combined_mesh.uv=coil_parts(1).coil_mesh.uv;
combined_mesh.boundary=coil_parts(1).coil_mesh.boundary;
combined_mesh.mesh_part_vertex_ind=ones(1,size(coil_parts(1).coil_mesh.vertices,2));
for part_ind=2:numel(coil_parts)
combined_mesh.faces=[combined_mesh.faces coil_parts(part_ind).coil_mesh.faces+size(combined_mesh.vertices,2)];
combined_mesh.n=[combined_mesh.n coil_parts(part_ind).coil_mesh.n];
combined_mesh.uv=[combined_mesh.uv coil_parts(part_ind).coil_mesh.uv];
combined_mesh.mesh_part_vertex_ind=[combined_mesh.mesh_part_vertex_ind ones(1,size(coil_parts(part_ind).coil_mesh.vertices,2))*part_ind];
for boundary_ind=1:numel(coil_parts(part_ind).coil_mesh.boundary  )
combined_mesh.boundary={combined_mesh.boundary{:} coil_parts(part_ind).coil_mesh.boundary{boundary_ind}+size(combined_mesh.vertices,2)};
end   
combined_mesh.vertices=[combined_mesh.vertices coil_parts(part_ind).coil_mesh.vertices];
end


%assign the stream function to the different mesh parts
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).stream_function=stream_function(coil_parts(part_ind).coil_mesh.unique_vert_inds);
end



end
