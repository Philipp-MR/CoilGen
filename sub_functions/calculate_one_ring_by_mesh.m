function coil_parts=calculate_one_ring_by_mesh(coil_parts,input)
%@Philipp Amrein 2021



for part_ind=1:numel(coil_parts)

coil_parts(part_ind).coil_mesh.faces=coil_parts(part_ind).coil_mesh.faces';
    
if ~input.temp_evalution.use_preoptimization_temp

num_nodes=size(coil_parts(part_ind).coil_mesh.vertices,2);
node_triangles = vertexAttachments(triangulation(coil_parts(part_ind).coil_mesh.faces,coil_parts(part_ind).coil_mesh.vertices'));
node_triangles_corners=cellfun(@(x) coil_parts(part_ind).coil_mesh.faces(x,:),node_triangles,'UniformOutput',0);
one_ring_list=cell(num_nodes,1);
one_ring_list=cell(num_nodes,1);%the triangles in which point hhh is part of

for node_ind=1:num_nodes
single_cell=num2cell(node_triangles_corners{node_ind},2);
one_ring_list{node_ind}=cell2mat(cellfun(@(x) x(x~=node_ind),single_cell,'UniformOutput',0))';
end



%make sure that the current orientation is uniform for all elements  (old)
for node_ind=1:num_nodes 
for face_ind=1:size(one_ring_list{node_ind},2)
point_aa=[coil_parts(part_ind).coil_mesh.vertices(:,one_ring_list{node_ind}(1,face_ind))];
point_bb=[coil_parts(part_ind).coil_mesh.vertices(:,one_ring_list{node_ind}(2,face_ind))];
point_cc=[coil_parts(part_ind).coil_mesh.vertices(:,node_ind)];   
cross_vec=cross(point_bb-point_aa,point_aa-point_cc);
if sign(dot(coil_parts(part_ind).coil_mesh.n(:,node_ind)',cross_vec))>0
one_ring_list{node_ind}(:,face_ind)=flipud(one_ring_list{node_ind}(:,face_ind));
end
end
end


%order the current element in circular arragnment (old)
node_triangle_mat=false(num_nodes,size(coil_parts(part_ind).coil_mesh.faces,1));
% for node_ind=1:num_nodes
% %find the starting point in case for the boundary
% open_start=find(~arrayfun(@(x) any(one_ring_list{node_ind}(2,:)==x),one_ring_list{node_ind}(1,:)));
% if isempty(open_start) 
% open_start=1; 
% end
% cell_order=[open_start];
% next_el=one_ring_list{node_ind}(2,cell_order(end));
% while numel(cell_order)~=numel(one_ring_list{node_ind}(1,:))
% cell_order=[cell_order find(one_ring_list{node_ind}(1,:)==next_el)];
% next_el=one_ring_list{node_ind}(2,cell_order(end));
% end
% one_ring_list{node_ind}=one_ring_list{node_ind}(:,cell_order);
% node_triangles{node_ind}=node_triangles{node_ind}(cell_order); % also update the order for the one_ring_list
% node_triangle_mat(node_ind,node_triangles{node_ind})=true;
% end

coil_parts(part_ind).one_ring_list=one_ring_list;
coil_parts(part_ind).node_triangles=node_triangles;
coil_parts(part_ind).node_triangle_mat=node_triangle_mat;
coil_parts(part_ind).coil_mesh.faces=coil_parts(part_ind).coil_mesh.faces';

else

coil_parts(part_ind).one_ring_list=input.temp.coil_parts(part_ind).one_ring_list;
coil_parts(part_ind).node_triangles=input.temp.coil_parts(part_ind).node_triangles;
coil_parts(part_ind).node_triangle_mat=input.temp.coil_parts(part_ind).node_triangle_mat;
coil_parts(part_ind).coil_mesh.faces=coil_parts(part_ind).coil_mesh.faces';

end

end

% % %@Philipp Amrein 2019  OLD VERSION
% % num_nodes=size(coil_mesh.vertices,2);
% % %get all the vectors if the opposite edges for every node on the mesh
% % node_triangles=cell(num_nodes,1);%the triangles in which point hhh is part of
% % neighbour_point_list=cell(num_nodes,1);
% % for node_ind=1:num_nodes
% % node_triangles{node_ind}=find(any(coil_mesh.faces==node_ind,2));
% % for face_ind=1:numel(node_triangles{node_ind})
% % neighbour_point_list{node_ind}=[neighbour_point_list{node_ind} coil_mesh.faces(node_triangles{node_ind}(face_ind),:)];
% % end
% % neighbour_point_list{node_ind}=unique(neighbour_point_list{node_ind});
% % neighbour_point_list{node_ind}= neighbour_point_list{node_ind}(neighbour_point_list{node_ind}~=node_ind);
% % end
% % %find the right connectivity between the neighbours
% % neighbour_point_edge_list=cell(num_nodes,1);
% % for node_ind=1:num_nodes
% % for face_ind=1:numel(neighbour_point_list{node_ind})
% % neighbour_point_edge_list{node_ind}{face_ind}=intersect(neighbour_point_list{node_ind},neighbour_point_list{neighbour_point_list{node_ind}(face_ind)});
% % end
% % end
% %  %Find the opposing edges of triangles representated by a double node column
% % one_ring_list=cell(num_nodes,1);
% % for node_ind=1:num_nodes
% % num_edges=sum(cellfun(@numel,neighbour_point_edge_list{node_ind}));
% % one_ring_list{node_ind}=zeros(2,num_edges); % the actual indices of point indices that will define the current elements around the point of investigation
% % kkk=1;
% % for face_ind=1:numel(neighbour_point_edge_list{node_ind})
% % for fff=1:numel(neighbour_point_edge_list{node_ind}{face_ind})
% % one_ring_list{node_ind}(:,kkk)=[neighbour_point_list{node_ind}(face_ind) neighbour_point_edge_list{node_ind}{face_ind}(fff)];
% % kkk=kkk+1;
% % end
% % end
% % end
% % % delete double counted current element
% % for node_ind=1:num_nodes
% % kkk=0;
% % for fff=1:size(one_ring_list{node_ind},2)
% % if any(find((one_ring_list{node_ind}(1,:)==one_ring_list{node_ind}(2,fff-kkk))...
% % .*(one_ring_list{node_ind}(2,:)==one_ring_list{node_ind}(1,fff-kkk))))
% % one_ring_list{node_ind}(:,find((one_ring_list{node_ind}(1,:)==one_ring_list{node_ind}(2,fff-kkk))...
% % .*(one_ring_list{node_ind}(2,:)==one_ring_list{node_ind}(1,fff-kkk))))=[]; %#ok
% % kkk=kkk+1;
% % end
% % end
% % end
% % 
% % 
% % 
% % %order the current element in circular arragnment
% % for node_ind=1:num_nodes
% % %find the starting point in case for the boundary
% % open_start_elements=one_ring_list{node_ind}(arrayfun(@(x) sum(one_ring_list{node_ind}(:)==x,'all'),[one_ring_list{node_ind}(:)])==1); %check if some nodes are part of an open end
% % if isempty(open_start_elements) 
% % start_edge=1;  %edges form closed loop
% % else
% % start_edge=find(arrayfun(@(x) any(one_ring_list{node_ind}(:,x)==open_start_elements(1)),1:size(one_ring_list{node_ind},2)));
% % %Flip the start edge if it has the wrong orientation
% % if one_ring_list{node_ind}(2,start_edge)==open_start_elements(1)
% % one_ring_list{node_ind}(:,start_edge)=flipud(one_ring_list{node_ind}(:,start_edge));
% % end
% % end
% % cell_order=[start_edge];
% % next_el=one_ring_list{node_ind}(2,start_edge);
% % while numel(cell_order)~=numel(one_ring_list{node_ind}(1,:))
% % next_pos_edge_index=setdiff(find(arrayfun(@(x) any(one_ring_list{node_ind}(:,x)==next_el),1:size(one_ring_list{node_ind},2))),cell_order(end)); %find the edge indexs wich carry the next element but which is not the last
% % %Flip the next edge if it has the wrong orientation
% % if one_ring_list{node_ind}(2,next_pos_edge_index)==next_el
% % one_ring_list{node_ind}(:,next_pos_edge_index)=flipud(one_ring_list{node_ind}(:,next_pos_edge_index));
% % end
% % cell_order=[cell_order next_pos_edge_index];
% % next_el=one_ring_list{node_ind}(2,cell_order(end));
% % end
% % one_ring_list{node_ind}=one_ring_list{node_ind}(:,cell_order);
% % node_triangles{node_ind}=node_triangles{node_ind}(cell_order); % also update the order for the one_ring_list
% % end
% % 
% % 
% % %make sure that the current orientation is uniform for all elements, check
% % %check the chirality with mesh normals
% % for node_ind=1:num_nodes 
% % vec_a=coil_mesh.vertices(:,one_ring_list{node_ind}(1,1))-coil_mesh.vertices(:,node_ind);
% % vec_b=coil_mesh.vertices(:,one_ring_list{node_ind}(1,2))-coil_mesh.vertices(:,one_ring_list{node_ind}(1,1));
% % cross_vec=cross(vec_a,vec_b);
% % if dot(cross_vec,coil_mesh.n(:,node_ind))<0 
% % one_ring_list{node_ind}=fliplr(one_ring_list{node_ind});
% % one_ring_list{node_ind}=flipud(one_ring_list{node_ind});
% % node_triangles{node_ind}=flipud(node_triangles{node_ind});
% % end
% % end




end