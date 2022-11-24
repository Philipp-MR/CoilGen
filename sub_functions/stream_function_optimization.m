function [coil_parts,combined_mesh,b_field_opt_sf] = stream_function_optimization(coil_parts,target_field,input)
tikonov_reg_factor=input.tikonov_reg_factor;
%combine the matrices from the different mesh parts
sensitivity_matrix=[];
gradient_sensitivity_matrix=[];
resistance_matrix=[];
current_density_mat=[];

if ~input.temp_evalution.use_optimized_temp

%Concatinate the matrices over the different mesh parts
for part_ind=1:numel(coil_parts)
if part_ind==1
current_density_mat=coil_parts(part_ind).current_density_mat;
else
current_density_mat=cat(3 ...
    ,blkdiag(squeeze(current_density_mat(:,:,1)),squeeze(coil_parts(part_ind).current_density_mat(:,:,1)))...
    ,blkdiag(squeeze(current_density_mat(:,:,2)),squeeze(coil_parts(part_ind).current_density_mat(:,:,2)))...
    ,blkdiag(squeeze(current_density_mat(:,:,3)),squeeze(coil_parts(part_ind).current_density_mat(:,:,3))));
end
sensitivity_matrix=cat(3,sensitivity_matrix,coil_parts(part_ind).sensitivity_matrix);
gradient_sensitivity_matrix=cat(3,gradient_sensitivity_matrix,coil_parts(part_ind).sensitivity_matrix);
resistance_matrix=blkdiag(resistance_matrix,coil_parts(part_ind).resistance_matrix);
end

%generate a combined mesh container
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
combined_mesh.bounding_box=[min(combined_mesh.vertices(1,:)) max(combined_mesh.vertices(1,:)); min(combined_mesh.vertices(2,:)) max(combined_mesh.vertices(2,:)); min(combined_mesh.vertices(3,:)) max(combined_mesh.vertices(3,:))];
set_zero_flag=false; %flag to force the potential on the boundary nodes to zero
%Redece target field only to z component
sensitivity_matrix_single=[ squeeze(sensitivity_matrix(3,:,:))];
% if strcmp(input.target_field_definition_file,'x')
% gradient_sensitivity_matrix=[ squeeze(gradient_sensitivity_matrix(1,:,:))];
% elseif strcmp(input.target_field_definition_file,'y')
% gradient_sensitivity_matrix=[ squeeze(gradient_sensitivity_matrix(2,:,:))];
% else
% gradient_sensitivity_matrix=[ squeeze(gradient_sensitivity_matrix(3,:,:))];
% end
target_field_single=[ squeeze(target_field.b(3,:))];
%Reduce the Resistance matrix for boundary nodes
[reduced_res_matrix,~,~]= reduce_matrices_for_boundary_nodes(resistance_matrix,combined_mesh,set_zero_flag);
%Reduce the sensitivity matrix for boundary nodes
[reduced_sensitivity_matrix,boundary_nodes,is_not_boundary_node]= reduce_matrices_for_boundary_nodes(sensitivity_matrix_single,combined_mesh,set_zero_flag);
[reduced_gradient_sensitivity_matrix_x,~,~]= reduce_matrices_for_boundary_nodes(squeeze(gradient_sensitivity_matrix(1,:,:)),combined_mesh,set_zero_flag);
[reduced_gradient_sensitivity_matrix_y,~,~]= reduce_matrices_for_boundary_nodes(squeeze(gradient_sensitivity_matrix(2,:,:)),combined_mesh,set_zero_flag);
[reduced_gradient_sensitivity_matrix_z,~,~]= reduce_matrices_for_boundary_nodes(squeeze(gradient_sensitivity_matrix(3,:,:)),combined_mesh,set_zero_flag);
%Reduce the current density matrix for boundary nodes
[red_current_density_mat_u,~,~]= reduce_matrices_for_boundary_nodes(squeeze(current_density_mat(:,:,1))',combined_mesh,set_zero_flag);
[red_current_density_mat_v,~,~]= reduce_matrices_for_boundary_nodes(squeeze(current_density_mat(:,:,2))',combined_mesh,set_zero_flag);
[red_current_density_mat_w,~,~]= reduce_matrices_for_boundary_nodes(squeeze(current_density_mat(:,:,3))',combined_mesh,set_zero_flag);
%Scale the tikonov regularization factor with the number of target points and mesh vertices
tikonov_reg_factor=tikonov_reg_factor*size(reduced_sensitivity_matrix,1)/size(reduced_sensitivity_matrix,2);
if strcmp(input.sf_opt_method,'tikkonov')
%Calculate the stream function by the tikonov optimization approach
tik_reg_mat=tikonov_reg_factor*reduced_res_matrix;
reduced_sf=pinv(reduced_sensitivity_matrix'*reduced_sensitivity_matrix+tik_reg_mat'*tik_reg_mat)*reduced_sensitivity_matrix'*target_field_single';
else
%for initialization, calculate the tikkonov solution; then do a iteriative optimization
tik_reg_mat=tikonov_reg_factor*reduced_res_matrix;
% target_gradient_x=[ squeeze(target_field.target_gradient_dbdxyz(1,:))];
% target_gradient_y=[ squeeze(target_field.target_gradient_dbdxyz(2,:))];
% target_gradient_z=[ squeeze(target_field.target_gradient_dbdxyz(3,:))];
reduced_sf=pinv(reduced_sensitivity_matrix'*reduced_sensitivity_matrix+tik_reg_mat'*tik_reg_mat)*reduced_sensitivity_matrix'*target_field_single';
%reduced_sf=zeros(size(reduced_res_matrix,1),1);
%find the constrained solution
stream_func_max=max(reduced_sf)*2;
lb=ones(size(reduced_sf)).*(-1).*stream_func_max;
ub=ones(size(reduced_sf)).*stream_func_max;

%max_current=max(((red_current_density_mat_u*reduced_sf).^2+(red_current_density_mat_v*reduced_sf).^2+(red_current_density_mat_w*reduced_sf).^2).^(1/2));


cost_function = @(x) sum((reduced_sensitivity_matrix*x-target_field_single').^2)+tikonov_reg_factor*(x'*reduced_res_matrix*x);


% cost_function = @(x) sum((reduced_gradient_sensitivity_matrix_x*x-target_gradient_x').^2)+...
%                                     sum((reduced_gradient_sensitivity_matrix_y*x-target_gradient_y').^2)+...
%                                     sum((reduced_gradient_sensitivity_matrix_z*x-target_gradient_z').^2)+...
%     +(-1)*min((abs(reduced_gradient_sensitivity_matrix_x*x)))...
%     +tikonov_reg_factor*max(((red_current_density_mat_u*x).^2+(red_current_density_mat_v*x).^2+(red_current_density_mat_w*x).^2).^(1/2));


% cost_function = @(x) (-1)*min((abs(reduced_gradient_sensitivity_matrix_x*x)))...
%     +tikonov_reg_factor*max(((red_current_density_mat_u*x).^2+(red_current_density_mat_v*x).^2+(red_current_density_mat_w*x).^2).^(1/2));

%cost_function = @(x) (-1)*min((abs(reduced_gradient_sensitivity_matrix_x*x)))+sum((reduced_gradient_sensitivity_matrix_x*x-target_gradient_x').^2)+tikonov_reg_factor*(x'*reduced_res_matrix*x);
%cost_function = @(x) sum((reduced_gradient_sensitivity_matrix_x*x-target_gradient_x').^2)+tikonov_reg_factor*(x'*reduced_res_matrix*x);

%cost_function = @(x) max(reduced_gradient_sensitivity_matrix_x*x-target_gradient_x')+tikonov_reg_factor*(x'*reduced_res_matrix*x);

% cost_function = @(x) sum((reduced_gradient_sensitivity_matrix_x*x-target_gradient_x').^2)+...
%                                 sum((reduced_gradient_sensitivity_matrix_y*x-target_gradient_y').^2)+...
%                                 sum((reduced_gradient_sensitivity_matrix_z*x-target_gradient_z').^2)+...
%                                 tikonov_reg_factor*max(((red_current_density_mat_u*x).^2+(red_current_density_mat_v*x).^2+(red_current_density_mat_w*x).^2).^(1/2));

%options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options.MaxIterations=input.fmincon_parameter(1);
options.MaxFunctionEvaluations=input.fmincon_parameter(2);
options.OptimalityTolerance = input.fmincon_parameter(3);
options.ConstraintTolerance = input.fmincon_parameter(4);
options.StepTolerance = input.fmincon_parameter(5);
reduced_sf = fmincon(cost_function,...
                                                reduced_sf,...
                                                [], ...
                                                [], ...
                                                [], ...
                                                [], ...
                                                lb, ...
                                                ub, ...
                                                [], ...
                                                options);
end
%Reexpand the stream potential to the boundary nodes
opt_stream_func= reexpand_stream_function_for_boundary_nodes(reduced_sf,boundary_nodes,is_not_boundary_node,set_zero_flag);
combined_mesh.stream_function=opt_stream_func;
%Calculate the magnetic field generated by the optimized stream function
b_field_opt_sf=[squeeze(sensitivity_matrix(1,:,:))*opt_stream_func squeeze(sensitivity_matrix(2,:,:))*opt_stream_func squeeze(sensitivity_matrix(3,:,:))*opt_stream_func]';
%Seperate the optimized stream function again onto the different mesh parts
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).stream_function=opt_stream_func(combined_mesh.mesh_part_vertex_ind==part_ind);
% coil_parts(part_ind).current_density=[coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,combined_mesh.mesh_part_vertex_ind==part_ind,1); ...
%                                                         coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,combined_mesh.mesh_part_vertex_ind==part_ind,2); ...
%                                                         coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,combined_mesh.mesh_part_vertex_ind==part_ind,3)];
jx=coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,:,1);
jy=coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,:,2);
jz=coil_parts(part_ind).stream_function'*coil_parts(part_ind).current_density_mat(:,:,3);

coil_parts(part_ind).current_density=[jx; jy; jz];
end
% %Quadratic programming
% stream_func_max=3000;
% lb=ones(size(reduced_sf)).*(-1).*stream_func_max;
% ub=ones(size(reduced_sf)).*stream_func_max;
% H=reduced_sensitivity_matrix'*reduced_sensitivity_matrix;
% f=(-1).*reduced_sensitivity_matrix'*target_field_single';
% options.MaxIterations = 2.000000e+04;
% reduced_sf = quadprog(H,f,[],[],[],[],lb,ub,reduced_sf,options); %output of quadprog
% % %calculate the resulting current dentsity for the mesh faces

else

for part_ind=1:numel(coil_parts)
coil_parts(part_ind).current_density=input.temp.coil_parts(part_ind).current_density;
coil_parts(part_ind).stream_function=input.temp.coil_parts(part_ind).stream_function;
end

end

end
function [reduced_mat,boundary_nodes,is_not_boundary_node]= reduce_matrices_for_boundary_nodes(full_mat,coil_mesh,set_zero_flag)
%Reduce the sensitivtity matrix in order to limit the degrees of freedom on
%the boundary nodes and make sure that they have constant sf later for each
%boundary
num_nodes=size(coil_mesh.vertices,2);
dim_to_reduce=size(full_mat)==num_nodes;
num_boundaries = numel(coil_mesh.boundary);
num_nodes_per_boundary=arrayfun(@(x) numel(unique(coil_mesh.boundary{x})),1:num_boundaries);
% if num_boundaries==1
% num_nodes_per_boundary=arrayfun(@(x) numel(coil_mesh.boundary{x}),1:num_boundaries);
% else
% num_nodes_per_boundary=arrayfun(@(x) numel(coil_mesh.boundary{x}),1:num_boundaries)-1;
% end
is_not_boundary_node=setdiff(1:size(full_mat,2),cat(1,coil_mesh.boundary{:}));
boundary_nodes=arrayfun(@(x) unique(coil_mesh.boundary{x}),1:num_boundaries,'UniformOutput',false);
% reduced_mat_size=ones(1,ndims(full_mat)).*num_nodes;
% reduced_mat_size(dim_to_reduce)=size(full_mat,2)-sum(num_nodes_per_boundary)+num_boundaries;
reduced_mat=full_mat;
%Generate the summed contribution as the first entries
%for dim_to_reduce_ind=find(dim_to_reduce)
%Shift the dimension which will be reduced in first place
%reduced_mat  = permute(reduced_mat, [dim_to_reduce_ind, setdiff(1:ndims(full_mat), dim_to_reduce_ind)]);
if sum(dim_to_reduce==1)~=0
for dim_to_reduce_ind=find(dim_to_reduce)
Index1    = cell(1, ndims(full_mat));
Index1(:) = {':'};
Index2=Index1;
for boundary_ind=1:num_boundaries
if set_zero_flag==true 
Index1{dim_to_reduce_ind} = boundary_nodes{boundary_ind}(1);
reduced_mat(Index1{:})=0; 
else
Index1{dim_to_reduce_ind} = boundary_nodes{boundary_ind}(1);
Index2{dim_to_reduce_ind} = boundary_nodes{boundary_ind}(:);
reduced_mat(Index1{:})=sum(reduced_mat(Index2{:}),dim_to_reduce_ind);
end
end
end
%Reange the matrix to its reduced form
boundary_nodes_first_inds=arrayfun(@(x) boundary_nodes{x}(1),1:numel(boundary_nodes));
for dim_to_reduce_ind=find(dim_to_reduce)
prev_reduced_mat=reduced_mat;
Index1    = cell(1, ndims(full_mat));
Index1(:) = {':'};
Index1{dim_to_reduce_ind}=1:num_boundaries; %the first entries are the subsitutes for the boundaries
Index2    = cell(1, ndims(full_mat));
Index2(:) = {':'};
Index2{dim_to_reduce_ind}=boundary_nodes_first_inds; %these are the inds of the first nodes of the boundaries
Index3    = cell(1, ndims(full_mat));
Index3(:) = {':'};
Index3{dim_to_reduce_ind}=(num_boundaries+1):(numel(is_not_boundary_node)+num_boundaries); % these will be the inds of the non-boundary nodes
Index4    = cell(1, ndims(full_mat));
Index4(:) = {':'};
Index4{dim_to_reduce_ind}=is_not_boundary_node; % the old non boundary nodes
Index5    = cell(1, ndims(full_mat));
Index5(:) = {':'};
Index5{dim_to_reduce_ind}=((num_nodes-(sum(num_nodes_per_boundary)-num_boundaries))+1):num_nodes; %these entries will be deleted
reduced_mat(Index1{:})=prev_reduced_mat(Index2{:});
reduced_mat(Index3{:})=prev_reduced_mat(Index4{:});
reduced_mat(Index5{:})=[];
end
end
end
function sf= reexpand_stream_function_for_boundary_nodes(reduced_sf,boundary_nodes,is_not_boundary_node,set_zero_flag)
%reexpand the stream function to all nodes including the boundary nodes
%THE NODES OF THE FIRST BOUNDARY HAVE A POTENTIAL OF ZERO
sf=zeros(numel(is_not_boundary_node)+sum(arrayfun(@(x) numel(boundary_nodes{x}),1:numel(boundary_nodes))),1);
%assign the sf vals for the boundary nodes
for boundary_ind=1:numel(boundary_nodes)
if set_zero_flag==1
sf(boundary_nodes{boundary_ind})=0;
else
sf(boundary_nodes{boundary_ind})=reduced_sf(boundary_ind);
end
end
%asign the rest
sf(is_not_boundary_node)=reduced_sf(numel(boundary_nodes)+1:end);
end
% % % OLD  function [reduced_sens_matrix,boundary_nodes,is_not_boundary_node]= reduce_matrices_for_boundary_nodes(sens_matrix,coil_mesh)
% % % %Reduce the sensitivtity matrix in order to limit the degrees of freedom on
% % % %the boundary nodes and make sure that they have constant sf later for each
% % % %boundary
% % % num_boundaries = numel(coil_mesh.boundary);
% % % num_nodes_per_boundary=arrayfun(@(x) numel(coil_mesh.boundary{x}),1:num_boundaries)-1;
% % % reduced_sens_matrix=zeros(size(sens_matrix,1),size(sens_matrix,2)-sum(num_nodes_per_boundary)+num_boundaries);
% % % is_not_boundary_node=setdiff(1:size(sens_matrix,2),cat(1,coil_mesh.boundary{:}));
% % % boundary_nodes=arrayfun(@(x) unique(coil_mesh.boundary{x}),1:num_boundaries,'UniformOutput',false);
% % % 
% % % %Generate the summed contribution as the first entries
% % % for boundary_ind=1:num_boundaries
% % % reduced_sens_matrix(:,boundary_ind)=sum(sens_matrix(:,boundary_nodes{boundary_ind}),2);
% % % end
% % % %add the rest non-boundary elements
% % % reduced_sens_matrix(:,num_boundaries+1:end)=sens_matrix(:,is_not_boundary_node);
% % % 
% % % end
% % % 
% % % 
% % % 
% % % function sf= reexpand_stream_function_for_boundary_nodes(reduced_sf,boundary_nodes,is_not_boundary_node)
% % % %reexpand the stream function to all nodes including the boundary nodes
% % % sf=zeros(numel(is_not_boundary_node)+sum(arrayfun(@(x) numel(boundary_nodes{x}),1:numel(boundary_nodes))),1);
% % % 
% % % %assign the sf vals for the boundary nodes
% % % for boundary_ind=1:numel(boundary_nodes)
% % %     sf(boundary_nodes{boundary_ind})=reduced_sf(boundary_ind);
% % % end
% % % %asign the rest
% % % sf(is_not_boundary_node)=reduced_sf(numel(boundary_nodes)+1:end);
% % % 
% % % end

% % %calculate the resulting current dentsity for the mesh faces (OLD)
% for part_ind=1:numel(coil_parts)
% pot_diffs1=repmat(coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(3,:))-coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(1,:)),[1 3])';
% pot_diffs2=repmat(coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(2,:))-coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(1,:)),[1 3])';
% pot_diffs3=repmat(coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(3,:))-coil_parts(part_ind).stream_function(coil_parts(part_ind).coil_mesh.faces(2,:)),[1 3])';
% edge_vecs1=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(3,:))-coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(1,:));
% edge_vecs2=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(2,:))-coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(1,:));
% edge_vecs3=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(3,:))-coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(2,:));
% %calculate the resulting current 
% coil_parts(part_ind).current_density=edge_vecs1.*pot_diffs1+edge_vecs2.*pot_diffs2+edge_vecs3.*pot_diffs3;
% end
% % %calculate the resulting current dentsity for the mesh nodes
% for part_ind=1:numel(coil_parts)
% coil_parts(part_ind).vertex_current_density=zeros(3,size(coil_parts(part_ind).coil_mesh.vertices,2));
% for node_ind=1:size(coil_parts(part_ind).coil_mesh.vertices,2)
% neigbour_nodes=coil_parts(part_ind).one_ring_list{3}(1,:);
% pot_diffs=coil_parts(part_ind).stream_function(neigbour_nodes)'-coil_parts(part_ind).stream_function(node_ind)';
% vec_diffs=coil_parts(part_ind).coil_mesh.vertices(:,neigbour_nodes)-coil_parts(part_ind).coil_mesh.vertices(:,node_ind);
% %calculate the resulting current 
% coil_parts(part_ind).vertex_current_density(:,node_ind)=sum(vec_diffs.*repmat(pot_diffs,[3 1]),2);
% end
% end