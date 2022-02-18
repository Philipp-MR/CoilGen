function boundary_criteria_matrix=calculate_boundary_criteria_matrix(coil_mesh,basis_elements)

% calculate the socalled boundary_criteria_matrix which enforces a constant
% stream function on the boundary nodes. This will forbide currents to
% leave the current carrying surface
% num_nodes=size(parameterized_mesh.v,1);
% boundary_criteria_matrix=zeros(num_nodes,num_nodes);
% 
% for boundary_ind=1:numel(parameterized_mesh.loops)
%     
%     
%     
% boundary_nodes= parameterized_mesh.loops{boundary_ind}';
% is_boundary=ismember(1:num_nodes,boundary_nodes);
% 
% 
% 
% end

% num_nodes=size(parameterized_mesh.v,1);
% 
% %mark the nodes that are boundary nodes
% boundary_nodes=[];
% for boundary_ind=1:numel(parameterized_mesh.loops)
%     boundary_nodes=[boundary_nodes parameterized_mesh.loops{boundary_ind}'];
% end
% is_boundary=ismember(1:num_nodes,boundary_nodes);
% %define the criteria that currents do not leave surface
% boundary_criteria=[];
% for boundary_ind=1:numel(parameterized_mesh.loops)
% boundary_criteria=[boundary_criteria; [parameterized_mesh.loops{boundary_ind}(2:end) parameterized_mesh.loops{boundary_ind}(1:end-1)]];
% end
% 
% for hhhh=1:size(boundary_criteria,1)
% boundary_criteria_matrix(hhhh,boundary_criteria(hhhh,1))=1;
% boundary_criteria_matrix(hhhh,boundary_criteria(hhhh,2))=-1;
% end




% calculate the socalled boundary_criteria_matrix which enforces a constant
% stream function on the boundary nodes. This will forbide currents to
% leave the current carrying surface
num_nodes=size(coil_mesh.vertices,2);
boundary_criteria_matrix=zeros(num_nodes,num_nodes);

%mark the nodes that are boundary nodes
boundary_nodes=[];
for boundary_ind=1:numel(coil_mesh.boundary)
    boundary_nodes=[boundary_nodes coil_mesh.boundary{boundary_ind}'];
end
is_boundary=ismember(1:num_nodes,boundary_nodes);
%define the criteria that currents do not leave surface
boundary_criteria=[];
for boundary_ind=1:numel(coil_mesh.boundary)
boundary_criteria=[boundary_criteria; [coil_mesh.boundary{boundary_ind}(2:end) coil_mesh.boundary{boundary_ind}(1:end-1)]];
end

for hhhh=1:size(boundary_criteria,1)
boundary_criteria_matrix(hhhh,boundary_criteria(hhhh,1))=1;
boundary_criteria_matrix(hhhh,boundary_criteria(hhhh,2))=-1;
end




end