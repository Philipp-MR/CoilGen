function coil_parts=calculate_basis_functions(coil_parts)
% create the basis funtion container which represents the current density 

%Initialize the outputs
coil_parts(numel(coil_parts)).is_real_triangle_mat=[];
coil_parts(numel(coil_parts)).triangle_corner_coord_mat=[];
coil_parts(numel(coil_parts)).current_mat=[];
coil_parts(numel(coil_parts)).area_mat=[];
coil_parts(numel(coil_parts)).face_normal_mat=[];
coil_parts(numel(coil_parts)).basis_elements=[];

for part_ind=1:numel(coil_parts)


num_nodes=size(coil_parts(part_ind).coil_mesh.vertices,2);

%create the container for the basis function
coil_parts(part_ind).basis_elements(num_nodes).stream_function_potential=[];
coil_parts(part_ind).basis_elements(num_nodes).triangles=[];
num_triangles_per_node=[cellfun(@(x) size(x,2),coil_parts(part_ind).one_ring_list)]';
for node_ind=1:num_nodes 
node_point=coil_parts(part_ind).coil_mesh.vertices(:,node_ind); 
%Assign also the triangle by their indices of the mesh faces
coil_parts(part_ind).basis_elements(node_ind).triangles=coil_parts(part_ind).node_triangles{node_ind};
%assign stream function potential to this basis element
coil_parts(part_ind).basis_elements(node_ind).stream_function_potential=0;
for tri_ind=1:num_triangles_per_node(node_ind)
point_b=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).one_ring_list{node_ind}(1,tri_ind));
point_c=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).one_ring_list{node_ind}(2,tri_ind));
coil_parts(part_ind).basis_elements(node_ind).one_ring=coil_parts(part_ind).one_ring_list{node_ind}';

%calculate the area of the triangle
coil_parts(part_ind).basis_elements(node_ind).area(tri_ind)=norm(cross(point_c-node_point,point_b-node_point))/2;


% %Distance node to current element:  d = norm(cross(v1-v2,pt-v2)) / norm(v1-v2);         
% dist_node_current_element=norm(cross(point_b-point_c,node_point-point_c))/norm(point_b-point_c);
%face normal of the triangle
coil_parts(part_ind).basis_elements(node_ind).face_normal(tri_ind,:)=cross(point_c-node_point,point_c-point_b)./norm(cross(point_c-node_point,point_c-point_b));
%Corner points ABC of the triangle
coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,:,:)= [node_point,point_b,point_c]; %1st ind: num triangle, 2nd ind: coords, 3th ind: corner ind
%calc the tangential current density of the triangle
coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,:)=(point_c-point_b)./(2*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind));
end
end


%Create outputs in matrix form
highst_triangle_count_per_node=max(arrayfun(@(x) numel(coil_parts(part_ind).basis_elements(x).area),1:numel(coil_parts(part_ind).basis_elements)));
is_real_triangle_mat=false(num_nodes,highst_triangle_count_per_node);
triangle_corner_coord_mat=zeros(num_nodes,highst_triangle_count_per_node,3,3);
face_normal_mat=zeros(num_nodes,highst_triangle_count_per_node,3);
current_mat=zeros(num_nodes,highst_triangle_count_per_node,3);
area_mat=zeros(num_nodes,highst_triangle_count_per_node);
for node_ind=1:num_nodes
is_real_triangle_mat(node_ind,1:numel(coil_parts(part_ind).basis_elements(node_ind).area))=true;
triangle_corner_coord_mat(node_ind,is_real_triangle_mat(node_ind,:),:,:)=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC;
current_mat(node_ind,is_real_triangle_mat(node_ind,:),:)=coil_parts(part_ind).basis_elements(node_ind).current;
area_mat(node_ind,is_real_triangle_mat(node_ind,:))=coil_parts(part_ind).basis_elements(node_ind).area;
face_normal_mat(node_ind,is_real_triangle_mat(node_ind,:),:)=coil_parts(part_ind).basis_elements(node_ind).face_normal;
end

coil_parts(part_ind).is_real_triangle_mat=is_real_triangle_mat;
coil_parts(part_ind).triangle_corner_coord_mat=triangle_corner_coord_mat;
coil_parts(part_ind).current_mat=current_mat;
coil_parts(part_ind).area_mat=area_mat;
coil_parts(part_ind).face_normal_mat=face_normal_mat;



% % % %unify the current and triangle orientation
% % % if dot(basis_elements(hhhh).face_normal(gggg,:),[1;0;0])<0
% % % %Corner points ABC of the triangle
% % % basis_elements(hhhh).triangle_points_ABC(gggg,:,:)= [node_point,point_b,point_c];
% % % %calc the tangential current density of the triangle 
% % % basis_elements(hhhh).current(gggg,:)=(point_b-point_c).*(1/norm(point_b-point_c)).*(1/dist_node_current_element);
% % % else 
% % % %Corner points ABC of the triangle
% % % basis_elements(hhhh).triangle_points_ABC(gggg,:,:)= [node_point,point_b,point_c];
% % % %calc the tangential current density of the triangle 
% % % basis_elements(hhhh).current(gggg,:)=(point_b-point_c).*(1/norm(point_b-point_c)).*(1/dist_node_current_element).*(-1);
% % % end

end

end