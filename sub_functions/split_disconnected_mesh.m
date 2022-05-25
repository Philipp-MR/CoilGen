function coil_parts = split_disconnected_mesh( coil_mesh_in )

%Split the mesh and the stream function if there are disconnected pieces
%such as shielded gradients

coil_mesh_in.faces=coil_mesh_in.faces';

vert_group = zeros(size(coil_mesh_in.faces,1),1,'uint32');
current_group = 0;
while any(vert_group==0)
current_group = current_group + 1;
next_face = find(vert_group==0,1,'first');
verts_to_sort = coil_mesh_in.faces(next_face,:);
while ~isempty(verts_to_sort)
availFaceInds = find(vert_group==0);
[availFaceSub, ~] = find(ismember(coil_mesh_in.faces(availFaceInds,:), verts_to_sort)); 
vert_group(availFaceInds(availFaceSub)) = current_group;
verts_to_sort = coil_mesh_in.faces(availFaceInds(availFaceSub),:);
end
end
num_vert_groups = current_group;
coil_parts(num_vert_groups).coil_mesh=[];
for current_group = 1:num_vert_groups
faces_of_group = coil_mesh_in.faces(vert_group==current_group,:);
[unqVertIds, ~, newVertIndices] = unique(faces_of_group);
coil_parts(current_group).coil_mesh.faces = reshape(newVertIndices,size(faces_of_group))';
coil_parts(current_group).coil_mesh.vertices = coil_mesh_in.vertices(:,unqVertIds);
coil_parts(current_group).coil_mesh.unique_vert_inds=unqVertIds;
end

end