function coil_parts=calculate_resistance_matrix(coil_parts,target_field,input)


gauss_order=input.gauss_order;
conductor_thickness=input.conductor_thickness;
specific_conductivity_copper=input.specific_conductivity_conductor;
material_factor=specific_conductivity_copper/conductor_thickness;

for part_ind=1:numel(coil_parts)


%Check wether the matrix was already calcuated in a previous iteration
hash_target = generate_DataHash(target_field);
hash_ccs = generate_DataHash(input.coil_mesh_file);

if isfile(strcat(input.output_directory,"\temp\hash_ccs_part"+num2str(part_ind)+".txt"))
fileID = fopen(strcat(input.output_directory,"\temp\hash_ccs_part"+num2str(part_ind)+".txt"),'r');
read_string=textscan(fileID,'%c');
hash_ccs_temp=read_string{:}';
fclose(fileID);
else
fileID = fopen(strcat(input.output_directory,"\temp\hash_ccs_part"+num2str(part_ind)+".txt"),'w');
fprintf(fileID,'%c',hash_ccs);
hash_ccs_temp=' ';
fclose(fileID);
end
if isfile(strcat(input.output_directory,"\temp\hash_target_part"+num2str(part_ind)+".txt"))
fileID = fopen(strcat(input.output_directory,"\temp\hash_target_part"+num2str(part_ind)+".txt"),'r');
read_string=textscan(fileID,'%c');
hash_target_temp=read_string{:}';
fclose(fileID);
else
fileID = fopen(strcat(input.output_directory,"\temp\hash_target_part"+num2str(part_ind)+".txt"),'w');
fprintf(fileID,'%c',hash_target);
hash_target_temp=' ';
fclose(fileID);
end

if ~(strcmp(hash_target,hash_target_temp)&strcmp(hash_ccs,hash_ccs_temp)) | ~isfile(strcat(input.output_directory,'\temp\resistance_mat_temp_part',num2str(part_ind),'.mat'))

num_nodes=numel(coil_parts(part_ind).basis_elements);
    
%Calculate the adjacency that marks mesh node neighbours with coordinate
node_adjacency_mat = false(num_nodes);
for tri_ind = 1:size(coil_parts(part_ind).coil_mesh.faces,2)
  node_adjacency_mat(coil_parts(part_ind).coil_mesh.faces(1,tri_ind), coil_parts(part_ind).coil_mesh.faces(2,tri_ind)) = true;
  node_adjacency_mat(coil_parts(part_ind).coil_mesh.faces(2,tri_ind), coil_parts(part_ind).coil_mesh.faces(3,tri_ind)) = true;
  node_adjacency_mat(coil_parts(part_ind).coil_mesh.faces(3,tri_ind), coil_parts(part_ind).coil_mesh.faces(1,tri_ind)) = true;
end
[vert1,vert2]=find(node_adjacency_mat);
%mesh_edges=unique([vert1 vert2],'rows');
mesh_edges=[vert1 vert2];
mesh_edges_non_unique=[1:num_nodes mesh_edges(:,1)'; 1:num_nodes mesh_edges(:,2)']';

node_adjacency_mat = node_adjacency_mat | node_adjacency_mat';
coil_parts(part_ind).node_adjacency_mat=node_adjacency_mat;

%Calculate the matrix of spatial distances for the neighboring vertices


% nodal_neighbor_distances=((repmat(coil_parts(part_ind).coil_mesh.vertices(1,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])+[repmat(coil_parts(part_ind).coil_mesh.vertices(1,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])]').^2+...
%                                         (repmat(coil_parts(part_ind).coil_mesh.vertices(2,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])+[repmat(coil_parts(part_ind).coil_mesh.vertices(2,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])]').^2+...
%                                         (repmat(coil_parts(part_ind).coil_mesh.vertices(3,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])+[repmat(coil_parts(part_ind).coil_mesh.vertices(3,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1])]').^2).^(1/2);

    
% % %%% Old non vectorized version
% Calculate the resistance matrix Rmn
resistance_matrix=zeros(num_nodes,num_nodes); %initialize the sensitivity matrix
%for each face select its nodes; 
%find within the basis elements the triangles for this node and face
% multiply the current and area values of the involved basis functions and
% add the values
for edge_ind=1:size(mesh_edges_non_unique,1)
node_ind1=mesh_edges_non_unique(edge_ind,1);
node_ind2=mesh_edges_non_unique(edge_ind,2);
overlapping_triangles=intersect(coil_parts(part_ind).basis_elements(node_ind1).triangles,coil_parts(part_ind).basis_elements(node_ind2).triangles);
resistance_sum=0;
if ~isempty(overlapping_triangles)
for overlapp_tri_ind=overlapping_triangles
first_node_triangle_positon=coil_parts(part_ind).basis_elements(node_ind1).triangles==overlapp_tri_ind;
second_node_triangle_positon=coil_parts(part_ind).basis_elements(node_ind2).triangles==overlapp_tri_ind;
triangle_area=coil_parts(part_ind).basis_elements(node_ind1).area(first_node_triangle_positon);
%calcualte for node pairs the mutual resistance
%calcuate currents from the perspective for each corners
primary_current=coil_parts(part_ind).basis_elements(node_ind1).current(first_node_triangle_positon,:);
secondary_current=coil_parts(part_ind).basis_elements(node_ind2).current(second_node_triangle_positon,:);
%resistance_sum=resistance_sum+abs(dot(primary_current,secondary_current)*(triangle_area)^2);
resistance_sum=resistance_sum+dot(primary_current,secondary_current)*(triangle_area)^2;
end
resistance_matrix(node_ind1,node_ind2)=resistance_sum;
end
end
%resistance_matrix=resistance_matrix.*material_factor;
resistance_matrix=resistance_matrix+resistance_matrix';
resistance_matrix=resistance_matrix.*material_factor;

coil_parts(part_ind).resistance_matrix=resistance_matrix;

save(strcat(input.output_directory,'\temp\resistance_mat_temp_part',num2str(part_ind),'.mat'),'resistance_matrix');

else

load(strcat(input.output_directory,'\temp\resistance_mat_temp_part',num2str(part_ind),'.mat'),'resistance_matrix');
coil_parts(part_ind).resistance_matrix=resistance_matrix;


end




% % Unfinished vectorized version
% % Calculate the resistance matrix Rmn
% num_nodes=numel(basis_elements);
% resistance_matrix=zeros(num_nodes,num_nodes); %initialize the sensitivity matrix
% material_factor=specific_conductivity_copper/conductor_thickness;
% %for each face select its nodes; 
% %find within the basis elements the triangles for this node and face
% % multiply the current and area values of the involved basis functions and
% % add the values
% for node_ind1=1:num_nodes
% overlapping_triangles=arrayfun(@(x) intersect(basis_elements(node_ind1).triangles,basis_elements(x).triangles,'stable'),1:num_nodes,'UniformOutput',0);
% overlapping_triangles_inds=find(~arrayfun(@(x) isempty(overlapping_triangles{x}),1:num_nodes));
% first_node_tri_inds=basis_elements(node_ind1).triangles;
%     for node_ind2=overlapping_triangles_inds
%     second_node_tri_inds=basis_elements(node_ind2).triangles;
%     first_node_overlap_positon=find(ismember(first_node_tri_inds,overlapping_triangles{node_ind2}));
%     second_node_overlap_positon=find(ismember(second_node_tri_inds,overlapping_triangles{node_ind2}));
% %     [~,first_node_overlap_positon]=ismember(first_node_tri_inds,overlapping_triangles{node_ind2});
% %     [~,second_node_overlap_positon]=ismember(second_node_tri_inds,overlapping_triangles{node_ind2});
% %     first_node_overlap_positon=first_node_overlap_positon(first_node_overlap_positon~=0);
% %     second_node_overlap_positon=second_node_overlap_positon(second_node_overlap_positon~=0);
%     overlapp_triangle_areas=basis_elements(node_ind1).area(first_node_overlap_positon);
%     %calcualte for node pairs the mutual resistance
%     %calcuate currents from the perspective for each corners
%     primary_current=basis_elements(node_ind1).current(first_node_overlap_positon,:);
%     secondary_current=basis_elements(node_ind2).current(second_node_overlap_positon,:);
%     %resistance_sum=resistance_sum+abs(dot(primary_current,secondary_current)*(triangle_area)^2);
%     resistance_sum=sum(dot(primary_current',secondary_current').*(overlapp_triangle_areas).^2);
%     resistance_matrix(node_ind1,node_ind2)=resistance_sum;
%     resistance_matrix(node_ind2,node_ind1)=resistance_sum;
%     end
% end
% resistance_matrix=resistance_matrix.*material_factor;








% % % %%% Old non vectorized version
% % Calculate the resistance matrix Rmn
% num_nodes=numel(basis_elements);
% resistance_matrix=zeros(num_nodes,num_nodes); %initialize the sensitivity matrix
% material_factor=specific_conductivity_copper/conductor_thickness;
% %for each face select its nodes; 
% %find within the basis elements the triangles for this node and face
% % multiply the current and area values of the involved basis functions and
% % add the values
% for node_ind1=1:num_nodes
% for node_ind2=1:num_nodes
% overlapping_triangles=intersect(basis_elements(node_ind1).triangles,basis_elements(node_ind2).triangles);
% resistance_sum=0;
% if ~isempty(overlapping_triangles)
% for overlapp_tri_ind=overlapping_triangles
% first_node_triangle_positon=basis_elements(node_ind1).triangles==overlapp_tri_ind;
% second_node_triangle_positon=basis_elements(node_ind2).triangles==overlapp_tri_ind;
% triangle_area=basis_elements(node_ind1).area(first_node_triangle_positon);
% %calcualte for node pairs the mutual resistance
% %calcuate currents from the perspective for each corners
% primary_current=basis_elements(node_ind1).current(first_node_triangle_positon,:);
% secondary_current=basis_elements(node_ind2).current(second_node_triangle_positon,:);
% %resistance_sum=resistance_sum+abs(dot(primary_current,secondary_current)*(triangle_area)^2);
% resistance_sum=resistance_sum+dot(primary_current,secondary_current)*(triangle_area)^2;
% end
% resistance_matrix(node_ind1,node_ind2)=resistance_sum;
% resistance_matrix(node_ind2,node_ind1)=resistance_sum;
% end
% end
% end
% %resistance_matrix=resistance_matrix.*material_factor;
% resistance_matrix=resistance_matrix.*material_factor;






end