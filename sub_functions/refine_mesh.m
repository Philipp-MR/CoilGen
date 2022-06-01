function coil_parts=refine_mesh(coil_parts,input)
%increase the resoltion of the mesh and interpolate the stream function


iteration_num_mesh_refinement=input.iteration_num_mesh_refinement;
sf_source_file=input.sf_source_file;

if strcmp(sf_source_file,'none')

for part_ind=1:numel(coil_parts)


subdivided_mesh=coil_parts(part_ind).coil_mesh;
subdivided_mesh.faces=subdivided_mesh.faces';
subdivided_mesh.vertices=subdivided_mesh.vertices';

for num_subdivision_sf=1:iteration_num_mesh_refinement

    

%upsample the stream function
%calc edge centers
coord_1_3=[   arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[1 3])],1)) ,1:size(subdivided_mesh.faces,1)); ...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[1 3])],2)) ,1:size(subdivided_mesh.faces,1));...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[1 3])],3)) ,1:size(subdivided_mesh.faces,1))];
coord_3_2=[   arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[3 2])],1)) ,1:size(subdivided_mesh.faces,1)); ...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[3 2])],2)) ,1:size(subdivided_mesh.faces,1));...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[3 2])],3)) ,1:size(subdivided_mesh.faces,1)) ];
coord_2_1=[   arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[2 1])],1)) ,1:size(subdivided_mesh.faces,1)); ...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[2 1])],2)) ,1:size(subdivided_mesh.faces,1)); ...
                        arrayfun(@(x) mean(subdivided_mesh.vertices([subdivided_mesh.faces(x,[2 1])],3)) ,1:size(subdivided_mesh.faces,1)) ];

% %calc tri center
% new_center_coords=[  ...
% arrayfun(@(x) mean(subdivided_mesh.uv([subdivided_mesh.faces(x,:)],1)) ,1:size(subdivided_mesh.faces,1)); ...
% arrayfun(@(x) mean(subdivided_mesh.uv([subdivided_mesh.faces(x,:)],2)) ,1:size(subdivided_mesh.faces,1)) ];
% %use the distances as inverse weigthing 
% center_old_node_distances=[vecnorm(new_center_coords-subdivided_mesh.uv([subdivided_mesh.faces(:,1)],:)'); ...
%                                                                 vecnorm(new_center_coords-subdivided_mesh.uv([subdivided_mesh.faces(:,2)],:)'); ...
%                                                                 vecnorm(new_center_coords-subdivided_mesh.uv([subdivided_mesh.faces(:,3)],:)')];
%                                                             
% center_to_old_node_weights=1./center_old_node_distances;
% center_to_old_node_weights=center_to_old_node_weights./sum(center_to_old_node_weights,1);
% %
% %caclulate the interpolated stream function values on the new nodes
% center_potential=sum(subdivided_sf.pot(subdivided_mesh.faces).*center_to_old_node_weights',2);
% %


%upsample the mesh
    all_coords=[subdivided_mesh.vertices; coord_1_3'; coord_3_2'; coord_2_1'];
    coord_ind_1=subdivided_mesh.faces(:,1);
    coord_ind_2=subdivided_mesh.faces(:,2);
    coord_ind_3=subdivided_mesh.faces(:,3);
    new_coord_inds_1_3=[1:size(coord_1_3,2)]'+size(subdivided_mesh.vertices,1);
    new_coord_inds_3_2=[1:size(coord_3_2,2)]'+(size(subdivided_mesh.vertices,1)+size(coord_1_3,2));
    new_coord_inds_2_1=[1:size(coord_2_1,2)]'+(size(subdivided_mesh.vertices,1)+size(coord_1_3,2)+size(coord_3_2,2));
    %build the new triangles
    new_tri_1=[coord_ind_1 new_coord_inds_1_3 new_coord_inds_2_1]; 
    new_tri_2=[new_coord_inds_1_3 coord_ind_3 new_coord_inds_3_2];
    new_tri_3=[new_coord_inds_3_2 coord_ind_2 new_coord_inds_2_1];
    new_tri_4=[new_coord_inds_1_3 new_coord_inds_3_2 new_coord_inds_2_1];
    new_tri=[new_tri_1; new_tri_2; new_tri_3; new_tri_4];
    %delete double counted nodes
    [all_coords_unique,~,ic] =unique(all_coords ,'rows','stable');
    ind_replace_list=[[1:size(ic,1)]' ic];
    %new_tri = changem(new_tri,ind_replace_list(:,2),ind_replace_list(:,1));
    [to_replace, by_what] = ismember(new_tri, ind_replace_list(:,1));
    new_tri(to_replace) = ind_replace_list(by_what(to_replace),2);
    subdivided_mesh.faces=new_tri;
    subdivided_mesh.vertices=all_coords_unique;
    
    
end

coil_parts(part_ind).coil_mesh.vertices=subdivided_mesh.vertices';
coil_parts(part_ind).coil_mesh.faces=subdivided_mesh.faces';
    
end
    
end

end
