function sensitivity_matrix=calculate_sensitivity_matrix_without_gauss(basis_elements,coil_mesh,target_points)


% calculate the sensitivity matrix cn

% calculate the weights and the test point for the gauss legendre
% integration on each triangle is done by gauss legendre integration
%[u_coord,v_coord,gauss_weight] = gauss_legendre_integration_points_triangle(gauss_order); % find uv-coords and weights for the gaus-legendre integration
%num_gauss_points=numel(gauss_weight); % the number of gaus points
biot_savart_coeff=10^(-7); % biot savart coefficient
num_nodes=numel(basis_elements);
num_target_points=size(target_points,2);
sensitivity_matrix=zeros(3,num_target_points,num_nodes); %initialize the sensitivity matrix
num_triangles_per_node=arrayfun(@(x) numel(basis_elements(x).area),1:numel(basis_elements));

x_target=target_points(1,:);
y_target=target_points(2,:);
z_target=target_points(3,:);

dCx=zeros(1,num_target_points);
dCy=zeros(1,num_target_points);
dCz=zeros(1,num_target_points);

for node_ind=1:num_nodes % iterate of number of nodes
for tri_ind=1:num_triangles_per_node(node_ind) % iterate of number of triangles of that node

    
    
node_point=coil_mesh.vertices(:,node_ind); 
point_b=coil_mesh.vertices(:,basis_elements(node_ind).one_ring(tri_ind,2));
point_c=coil_mesh.vertices(:,basis_elements(node_ind).one_ring(tri_ind,1));
    
    
%extract the coords of the triangle corners
x1=node_point(1);
y1=node_point(2);
z1=node_point(3);
x2=point_b(1);
y2=point_b(2);
z2=point_b(3);
x3=point_c(1);
y3=point_c(2);
z3=point_c(3);
    
    
% %extract the coords of the triangle corners
% x1=basis_elements(hhhh).triangle_points_ABC(gggg,1,1); %1st ind: num triangle, 2nd ind: coords, 3th ind: corner ind
% y1=basis_elements(hhhh).triangle_points_ABC(gggg,2,1);
% z1=basis_elements(hhhh).triangle_points_ABC(gggg,3,1);
% x2=basis_elements(hhhh).triangle_points_ABC(gggg,1,2);
% y2=basis_elements(hhhh).triangle_points_ABC(gggg,2,2);
% z2=basis_elements(hhhh).triangle_points_ABC(gggg,3,2);
% x3=basis_elements(hhhh).triangle_points_ABC(gggg,1,3);
% y3=basis_elements(hhhh).triangle_points_ABC(gggg,2,3);
% z3=basis_elements(hhhh).triangle_points_ABC(gggg,3,3);

x_tri=(x1+x2+x3)/3;
y_tri=(y1+y2+y3)/3;
z_tri=(z1+z2+z3)/3;

% calculate the distance norm for the biot savart law
distance_norm=((x_tri-x_target).^2+(y_tri-y_target).^2+(z_tri-z_target).^2).^(-3/2);

vx=basis_elements(node_ind).current(tri_ind,1);
vy=basis_elements(node_ind).current(tri_ind,2);
vz=basis_elements(node_ind).current(tri_ind,3);

%calcute the contribution to the sensitivity matrix and add them up
dCx=dCx+((-1).*vz.*(y_target-y_tri)+vy.*(z_target-z_tri)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind);
dCy=dCy+((-1).*vx.*(z_target-z_tri)+vz.*(x_target-x_tri)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind);
dCz=dCz+((-1).*vy.*(x_target-x_tri)+vx.*(y_target-y_tri)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind);

% %calcute the contribution to the sensitivity matrix and add them up
% dCx=dCx+((zgauss_in_uv-z_target).*basis_elements(node_ind).current(tri_ind,2)-(ygauss_in_uv-y_target).*basis_elements(node_ind).current(tri_ind,3)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCy=dCy+((xgauss_in_uv-x_target).*basis_elements(node_ind).current(tri_ind,3)-(zgauss_in_uv-z_target).*basis_elements(node_ind).current(tri_ind,1)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCz=dCz+((ygauss_in_uv-y_target).*basis_elements(node_ind).current(tri_ind,1)-(xgauss_in_uv-x_target).*basis_elements(node_ind).current(tri_ind,2)).*distance_norm.*2.*basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);

end
sensitivity_matrix(:,:,node_ind)=[dCx' dCy' dCz']'.*biot_savart_coeff;
end

end