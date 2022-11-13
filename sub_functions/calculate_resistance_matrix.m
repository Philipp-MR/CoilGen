function coil_parts=calculate_resistance_matrix(coil_parts,input)


gauss_order=input.gauss_order;
conductor_thickness=input.conductor_thickness;
specific_conductivity_copper=input.specific_conductivity_conductor;
material_factor=specific_conductivity_copper/conductor_thickness;

for part_ind=1:numel(coil_parts)

if ~input.temp_evalution.use_preoptimization_temp

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


else

coil_parts(part_ind).resistance_matrix=input.temp.coil_parts(part_ind).resistance_matrix;

end


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