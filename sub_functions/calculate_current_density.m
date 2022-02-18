function current_density_mat= calculate_current_density(coil_mesh,basis_elements,tri_membership_list)

%calcualte the current density by the mesh and stream function potential


%calculate the resulting current for the mesh faces
current_density_mat=zeros(size(coil_mesh.faces,1),numel(basis_elements),3);
for node_ind=1:numel(tri_membership_list)
    for tri_ind=1:numel(tri_membership_list{node_ind})
        current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,1)=current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,1)+basis_elements(node_ind).current(tri_ind,1);
        current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,2)=current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,2)+basis_elements(node_ind).current(tri_ind,2);
        current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,3)=current_density_mat(tri_membership_list{node_ind}(tri_ind),node_ind,3)+basis_elements(node_ind).current(tri_ind,3);
    end
end

% x_mat=squeeze(current_density_mat(:,:,1));
% y_mat=squeeze(current_density_mat(:,:,2));
% z_mat=squeeze(current_density_mat(:,:,3));

% %calculate the resulting current for the mesh faces
% opt_current=zeros(size(coil_mesh.faces,1),3);
% for hhhh=1:size(sensitivity_matrix,3)
%     for gggg=1:numel(tri_membership_list{hhhh})
%         opt_current(tri_membership_list{hhhh}(gggg),:)=opt_current(tri_membership_list{hhhh}(gggg),:)+basis_elements(hhhh).current(gggg,:)*opt_stream_func_it(hhhh);
%     end
% end


%%calculate the resulting current for the mesh faces
%current_density_mat=zeros(size(coil_mesh.ConnectivityList,1),numel(basis_elements));
% for hhhh=1:numel(current_density_mat)
%     for gggg=1:numel(tri_membership_list{hhhh})
%         current_density_mat(tri_membership_list{hhhh}(gggg),:)=current_density_mat(tri_membership_list{hhhh}(gggg),:)+basis_elements(hhhh).current(gggg,:)*opt_stream_func_it(hhhh);
%     end
% end




end

